

## What is this notebook about?
In this notebook I analyze RGI summary and read count files at multiple downsampling depths for *Acinetobacter baumannii* (and I did it with other pathogens too), clean and unify them, compute per-identifier resistant gene counts, merge with total reads to derive percent resistant reads, filter out predefined low-quality combinations, and map barcodes to hospital labels.
I apply a bootstrap routine to estimate medians and 95% confidence intervals of resistant proportions by downsampling level and identifier, aggregate counts for some relevant genes such as OXA-441 and also for all resistance genes, and produce summary tables of resistant fraction.
I visualize these results with line plots showing median trends over downsampling (GB), overlay bootstrap confidence ribbons and jittered replicate points, and distinguish hospitals using custom colors and shapes in separate plots for the single gene and the combined resistance profile.

## Understanding the notebook
Input: RGI summary and read count files for *Acinetobacter baumannii* at various downsampling levels, stored in rgi_concatenate and rgi_count_all_reads directories, respectively.
Output: See `plots` directory for the generated plots.

## Credits
The plot `downsampling_resistome_variation/plots/proportion_resistant_baumannii.png` and this script 
was generated together with [Diego Taquiri](https://github.com/diego-taquiri).

```{r load libraries}
library(tidyverse)
library(readr)
library(here)
library(janitor)
library(ggplot2)
library(scales)
library(boot)

```

Base functions:

```{r functions}
# Define a function to calculate variance, required by the boot function
calc_variance <- function(data, indices) {
  sampled_data <- data[indices]
  return(var(sampled_data))
}

# Function to calculate the median and its 95% CI using bootstrap
median_cl_boot <- function(x, conf = 0.95) {
  med <- median(x)
  
  # Check if all values are the same
  if (length(unique(x)) == 1) {
    return(data.frame(y = med, ymin = NA, ymax = NA))
  }
  
  # Bootstrap resampling
  boot_res <- boot(x, function(x, i) median(x[i]), R = 1000)
  
  # Check for error in bootstrapping
  if (all(boot_res$t == med)) {
    return(data.frame(y = med, ymin = NA, ymax = NA))
  }
  
  # Calculate confidence interval
  ci <- boot.ci(boot_res, type = "perc", conf = conf)$percent[4:5]
  
  # Return a data frame with the median and confidence interval
  data.frame(y = med, ymin = ci[1], ymax = ci[2])
}
```

Define themes and color schemes
```{r color scheme}
# Define the color scheme 
color_scheme <- c("Phenol"="#CAB2D6",
                  "Zymo Kit"="#1F78B4",
                  "Phenol + Frag"="#A5CEE3",
                  "Zymo Kit + Frag"="#6A3D9A")
color_scheme2 <- c("NT"="#F89C74", "80 ºC"="#F89897")

my_theme <- theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.y = element_text(size = 14),
        legend.position = "none", 
        plot.title.position = "plot",  
        plot.title = element_text(hjust = 0.5))
my_theme2 <- theme(
  axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
  axis.text.y = element_text(size = 12),
  axis.title.y = element_text(size = 14),
  axis.title.x = element_text(size = 14),
  plot.title.position = "plot",
  plot.title = element_text(hjust = 0.5, size=14),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.line = element_line(colour = "black", size = 0.5),
  legend.key = element_rect(fill = "white", colour = NA),
  strip.text = element_text(face = "bold", size = rel(1))
)

```

Read all the RGI from DS and fasta counts
```{r}
# Definir el vector de tamaños de downsampling y las rutas de archivos
downsampling_sizes <- c("03", "05", "1", "3", "6", "9")
file_path_reads <- "C:\\Users\\DAVID 21\\OneDrive\\Documentos\\Mirkoslab\\loui\\downsampling\\rgi_concatenate\\Acinetobacter_baumannii\\DS_%s_combined_rgi_summary.tsv"
file_path_counts <- "C:\\Users\\DAVID 21\\OneDrive\\Documentos\\Mirkoslab\\loui\\downsampling\\rgi_count_all_reads\\Acinetobacter_baumannii\\DS_%s_read_counts.txt"

# Crear una lista de data frames para reads y counts
DS_reads_list <- downsampling_sizes %>%
  map(~ {
    read_tsv(sprintf(file_path_reads, .x)) %>%
      clean_names() %>%
      mutate(downsamplingGB = ifelse(.x == "03", 0.3, 
                                     ifelse(.x == "05", 0.5, as.numeric(.x)))) %>%
      select(downsamplingGB, everything())
  })

DS_counts_list <- downsampling_sizes %>%
  map(~ {
    read_tsv(sprintf(file_path_counts, .x)) %>%
      clean_names() %>%
      mutate(downsamplingGB = ifelse(.x == "03", 0.3, 
                                     ifelse(.x == "05", 0.5, as.numeric(.x)))) %>%
      select(downsamplingGB, everything())
  })

# Asignar cada data frame a su respectiva variable en el entorno global
list2env(setNames(DS_reads_list, paste0("DS_", downsampling_sizes)), .GlobalEnv)
list2env(setNames(DS_counts_list, paste0("DS_", downsampling_sizes, "_count")), .GlobalEnv)
```

Unifiy all the tables of resistance
```{r}
all_downsamp <- bind_rows(DS_03, DS_05, DS_1, DS_3, DS_6, DS_9) %>%
  mutate(barcode_name=case_when(
    identifier%in%c("barcode01","barcode02","barcode03","barcode04")~"HMA AR070_1",
    identifier%in%c("barcode05","barcode06","barcode07")~"HMA AR070_2",
    identifier%in%c("barcode09","barcode10","barcode11","barcode12")~"HMA AR070_3",
    .default = identifier
  )) 
```

Get the count of resistance genes per downsampling
```{r}
# Creo mi lista con los nombres de los data frames
ds_names <- c("DS_03", "DS_05", "DS_1", "DS_3", "DS_6", "DS_9")

# Aplicar la operación de group_by, summarise, y ungroup
processed_dfs <- ds_names %>%
  map(~ get(.x) %>%   # Obtener el data frame usando su nombre
        group_by(downsamplingGB, identifier) %>% 
        summarise(resistant_count = n(), .groups = 'drop') %>% 
        ungroup())

# Asignar los data frames procesados a variables con sufijo "_summary" en el entorno global
list2env(setNames(processed_dfs, paste0(ds_names, "_table")), .GlobalEnv)

all_downsamp_table <- all_downsamp %>% 
  group_by(downsamplingGB,identifier) %>% 
  summarise(resistant_count=n()) %>% 
  ungroup()

```
Unify the counts of reads
```{r}
all_downsamp_count <- bind_rows(DS_03_count, DS_05_count, DS_1_count, DS_3_count, DS_6_count, DS_9_count)
```
Find the proportion of resistant reads and filter insufficient data in DSs
```{r}
all_DS_proportion <- all_downsamp_table %>% 
  left_join(all_downsamp_count,by=c("downsamplingGB","identifier")) %>% 
  mutate(porportion=resistant_count/read_count*100) %>% 
  separate(identifier,into = c("identifier","gene"),sep = "_") %>% 
  mutate(barcode_name=case_when(
    identifier%in%c("barcode01","barcode02","barcode03","barcode04")~"HMA AR070_1",
    identifier%in%c("barcode05","barcode06","barcode07")~"HMA AR070_2",
    identifier%in%c("barcode09","barcode10","barcode11","barcode12")~"HMA AR070_3",
    .default = identifier
  ))
all_DS_proportion

###FILTRAR PIPIPI
# Definir las combinaciones que queremos excluir
exclude_combinations <- tibble(
  downsamplingGB = c(9, 9, 9, 9, 6),
  identifier = c("barcode03", "barcode10", "barcode11", "barcode12", "barcode10")
)

# Aplicar anti_join para eliminar las combinaciones indeseadas
all_DS_proportion <- all_DS_proportion %>%
  anti_join(exclude_combinations, by = c("downsamplingGB", "identifier"))

# Mostrar el resultado filtrado
str(all_DS_proportion)

```

Make the graphs for each one

```{r}

oxa253_data <- all_DS_proportion %>% 
  filter(gene == "OXA-441") %>% group_by(downsamplingGB, identifier, barcode_name) %>% summarise(mean_proportion=mean(porportion)) #para que salga un punto por cada barcode, y no un punto por cada gen

p<-ggplot(oxa253_data, aes(x = downsamplingGB, y = mean_proportion)) +
  geom_line(aes(group = barcode_name, color = barcode_name), 
            stat = "summary", fun = median, size = 1)+
  labs(title = expression(paste("Resistant Proportion for Gene ", italic("adeN"))),
       x = "Downsampling Level (GB)",
       y = "Proportion") +
    scale_x_continuous(breaks = c(0.3,0.5,1,3,6,9),
                     labels = c(0.3,0.5,1,3,6,9))+
    # Add 95% confidence interval
  stat_summary(
    fun.data = median_cl_boot, 
    aes(group = barcode_name, fill = barcode_name), 
    geom = "ribbon", 
    alpha = 0.3
  )+  #geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color = barcode_name,shape = barcode_name),
              height = 0,width = 0.05, size = 2.5) +
  scale_color_manual(values = c(
    "HMA AR070_1" = "#E3211C",
    "HMA AR070_2" = "#0072B2",
    "HMA AR070_3" = "#F0E442"
  )) +
  scale_shape_manual(values = c(
    "HMA AR070_1" = 15,  # Filled triangle
    "HMA AR070_2" = 18,  # Filled square
    "HMA AR070_3" = 17  # Filled circle
  ))+
  scale_fill_manual(values = c(
    "HMA AR070_1" = "#E3211C",
    "HMA AR070_2" = "#0072B2",
    "HMA AR070_3" = "#F0E442"))+
  my_theme2 +
  guides(color = guide_legend(title = "Hospital"),
         shape = guide_legend(title = "Hospital"),
         group = guide_legend(title = "Hospital"),
         fill = guide_legend(title = "Hospital"))
```

Make the graph for all genes
```{r}

all_resistance <- all_DS_proportion %>% group_by(downsamplingGB, identifier, barcode_name) %>% summarise(mean_proportion=mean(porportion))

p<-ggplot(all_resistance, aes(x = downsamplingGB, y = mean_proportion)) +
  geom_line(aes(group = barcode_name, color = barcode_name), 
            stat = "summary", fun = median, size = 1)+
  labs(title = "Resistant Proportion for All Genes",
       x = "Downsampling Level (GB)",
       y = "Proportion") +
    scale_x_continuous(breaks = c(0.3,0.5,1,3,6,9),
                     labels = c(0.3,0.5,1,3,6,9))+
    # Add 95% confidence interval
  stat_summary(
    fun.data = median_cl_boot, 
    aes(group = barcode_name, fill = barcode_name), 
    geom = "ribbon", 
    alpha = 0.3
  )+  #geom_boxplot(outlier.shape = NA)+
  geom_jitter(aes(color = barcode_name,shape = barcode_name),
              height = 0,width = 0.05, size = 2.5) +
  scale_color_manual(values = c(
    "HMA AR070_1" = "#E3211C",
    "HMA AR070_2" = "#0072B2",
    "HMA AR070_3" = "#F0E442"
  )) +
  scale_shape_manual(values = c(
    "HMA AR070_1" = 15,  # Filled triangle
    "HMA AR070_2" = 18,  # Filled square
    "HMA AR070_3" = 17  # Filled circle
  ))+
  scale_fill_manual(values = c(
    "HMA AR070_1" = "#E3211C",
    "HMA AR070_2" = "#0072B2",
    "HMA AR070_3" = "#F0E442"))+
  my_theme2 +
  guides(color = guide_legend(title = "Hospital"),
         shape = guide_legend(title = "Hospital"),
         group = guide_legend(title = "Hospital"),
         fill = guide_legend(title = "Hospital"))
```




