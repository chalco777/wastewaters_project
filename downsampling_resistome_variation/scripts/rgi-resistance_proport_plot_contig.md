
## What is this notebook about?

In this notebook I estimate, across multiple downsampling (DS) depths, the **fraction of reads that map to AMR-labeled regions, specifically carbapenem-resistance gene (CRE) regions, and are themselves classified as carrying these ARGs**. Focusing on *Acinetobacter baumannii* (and extended to other pathogens), I:

* clean and unify the input summaries,
* compute per-AMR-region counts of reads classified as AMR,
* merge those with the total reads mapped to AMR regions to derive the percent of truly AMR reads,
*filter low-quality combinations, and
* map barcodes to their original hospital labels.

I then performed a bootstrap routine to obtain medians and 95% confidence intervals of the resistant proportion by downsampling level and identifier, aggregate results for relevant genes (e.g., OXA-441) as well as all the antibiotic resistance genes, and produce tables and figures of the resistant fraction. Plots show median trends vs. downsampling (GB), with bootstrap ribbons, jittered replicates, and hospital-specific colors/shapes—both for single genes and the combined AMR profile.

## Understanding the notebook
Input.
* RGI summaries run on **ONT reads** (not contigs), which label per-read ORFs/regions with AMR calls.
* Read-count files with the **number of reads that map to those AMR-labeled regions** (same identifiers), provided at several downsampling levels for A. baumannii and other pathogens. Files reside under rgi_concatenate (summaries) and rgi_count_all_reads (counts).

Output.
* See the `plots/` directory for generated figures and the summary tables of resistant fractions.

### Terminology (for clarity).
* AMR-mapped reads: reads that align to regions labeled as AMR by RGI.
* Truly AMR reads: among AMR-mapped reads, those that RGI classifies as AMR by themselves.
* Identifier: the per-region label (derived from the ONT read and ORF) used to join summaries and counts.


## Credits
The plot `downsampling_resistome_variation/plots/proportion_resistant_baumannii.png` and part of this script 
was generated together with [Diego Taquiri](https://github.com/diego-taquiri).

## Load libraries

```{r load libraries}
library(tidyverse)
library(readr)
library(here)
library(janitor)
library(ggplot2)
library(scales)
library(boot)

```

## Define Base functions:

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

## Define themes and color schemes
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

## Read the RGI results on AMR-mapped reads for each DS, and then the number of AMR-mapped reads for each Pathogen and DS

As expected, from these files I will get the proportion of trully carbapenem-resistant reads

```{r}
# Define the vector of downsampling sizes and file paths
downsampling_sizes <- c("03", "05", "1", "3", "6", "9")
file_path_reads <- "C:\\Users\\DAVID 21\\OneDrive\\Documentos\\Mirkoslab\\loui\\downsampling\\rgi_concatenate\\Acinetobacter_baumannii\\DS_%s_combined_rgi_summary.tsv"
file_path_counts <- "C:\\Users\\DAVID 21\\OneDrive\\Documentos\\Mirkoslab\\loui\\downsampling\\rgi_count_all_reads\\Acinetobacter_baumannii\\DS_%s_read_counts.txt"

# Create a list of data frames for reads and counts
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

# Assign each data frame to its respective variable in the global environment
list2env(setNames(DS_reads_list, paste0("DS_", downsampling_sizes)), .GlobalEnv)
list2env(setNames(DS_counts_list, paste0("DS_", downsampling_sizes, "_count")), .GlobalEnv)
```

## Unify all the tables of trully AMR genes (RGI output on reads)
```{r}
all_downsamp <- bind_rows(DS_03, DS_05, DS_1, DS_3, DS_6, DS_9) %>%
  mutate(barcode_name=case_when(
    identifier%in%c("barcode01","barcode02","barcode03","barcode04")~"HMA AR070_1",
    identifier%in%c("barcode05","barcode06","barcode07")~"HMA AR070_2",
    identifier%in%c("barcode09","barcode10","barcode11","barcode12")~"HMA AR070_3",
    .default = identifier
  )) 
```

## Get the count of resistant reads per gene from the RGI output
```{r}
# Create my list with the names of the data frames
ds_names <- c("DS_03", "DS_05", "DS_1", "DS_3", "DS_6", "DS_9")

# Apply the group_by, summarise, and ungroup operations
processed_dfs <- ds_names %>%
  map(~ get(.x) %>%   # Get the data frame using its name
        group_by(downsamplingGB, identifier) %>% 
        summarise(resistant_count = n(), .groups = 'drop') %>% 
        ungroup())

# Assign the processed data frames to variables with suffix "_summary" in the global environment
list2env(setNames(processed_dfs, paste0(ds_names, "_table")), .GlobalEnv)

all_downsamp_table <- all_downsamp %>% 
  group_by(downsamplingGB,identifier) %>% 
  summarise(resistant_count=n()) %>% 
  ungroup()

```
## Unify the counts of total mapped-reads to carbapenem resistant regions
```{r}
all_downsamp_count <- bind_rows(DS_03_count, DS_05_count, DS_1_count, DS_3_count, DS_6_count, DS_9_count)
```
## Find the proportion of trully resistant reads and filter insufficient data in downsamplings (DSs)
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

###FILTER PIPIPI
# Define the combinations we want to exclude
exclude_combinations <- tibble(
  downsamplingGB = c(9, 9, 9, 9, 6),
  identifier = c("barcode03", "barcode10", "barcode11", "barcode12", "barcode10")
)

# Apply anti_join to remove the unwanted combinations
all_DS_proportion <- all_DS_proportion %>%
  anti_join(exclude_combinations, by = c("downsamplingGB", "identifier"))

# Show the filtered result
str(all_DS_proportion)

```

## Make the graphs for each relevant CRE gene, changing the name of the gene in `filter(gene == "OXA-441")`

```{r}

oxa253_data <- all_DS_proportion %>% 
  filter(gene == "OXA-441") %>% group_by(downsamplingGB, identifier, barcode_name) %>% summarise(mean_proportion=mean(porportion)) #so that there is one point per barcode, and not one point per gene

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

## Make the graph for the proportion of resistant reads across all carbapenem-resistant genes
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





