#!/usr/bin/env Rscript
---
title: "pathogen_abundance"
author: "Adrián Chalco"
date: "2025-01-10"
output: html_document
---

What is this notebook about?

Here I assemble cross-run metadata (barcodes, place, date) and total read counts (MW, temporal, 17_set, seq2111, seq1212) to compute pathogen-level **relative abundance (%)** from Sylph outputs (with seq2111/seq1212 harmonized), keeping the legacy ConQuR path disabled. 
I exclude controls, reshape by species and place, and generate **temporal faceted plots** with median ± IQR as well as **library-size** trends (boxplots + median lines). I then overlay **crude clinical events** as date markers and quantify associations via two **correlation** figures (fraction vs total reads; assigned reads vs total, with linear fits and `ggpubr::stat_cor`), saving all panels as publication-ready PNGs.


## Preparar librerías
```{r}
library(tidyverse)
library(readr)
library(janitor)
library(ggplot2)
library(ggpubr)
library(scales)
library(boot)
library(RColorBrewer)
library(here)

```
## Preparar metadata para cada barcode de cada carpeta de secuenciamiento
```{r}
# Crear el data frame con los datos proporcionados
met_temp <- data.frame(
  `sample` = c("AR071_1", "AR071_2", "AR071_2", "AR071_2", "AR071_3",
               "AR072_1", "AR072_1", "AR072_1", "AR072_2", "AR072_2",
               "AR072_2", "AR072_3", "Mock Community", "Negative Control - Human DNA"),
  barcode = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
)

# Convertir la columna Barcode en el formato barcodeXX
met_temp$barcode <- sprintf("barcode%02d", met_temp$barcode)

# Mostrar el data frame resultante
print(met_temp)

# Crear el segundo data frame con los nuevos datos
additional_data <- data.frame(
  sample = c("AR071_1", "AR071_2", "AR071_3", 
             "AR072_1", "AR072_2", "AR072_3"),
  place = rep("HMA", 6),  # Columna 'HMA' con valor constante
  date = c("13/02/2024", "13/02/2024", "13/02/2024", 
           "20/02/2024", "20/02/2024", "20/02/2024")
)
# Unir ambos data frames por la columna 'sample_of_extraction'
temp_met <- merge(met_temp, additional_data, by = "sample", all.x = TRUE)

#Para L
# Crear el tercer data frame con los datos del 2do intento
el_met<- data.frame(
  sample = c("AR086", "AR093","AR093","AR093", "AR095", 
             "Mock Community", "Negative Control - Human DNA"),
  `date` = c("06/05/2024", "22/05/2024", "22/05/2024", "22/05/2024", "27/05/2024", NA, NA),
  place = c("INSN", "HMA", "HMA", "HMA", "INSN",NA, NA),
  `barcode` = c(24, 23, 22, 21, 20, 19, 18)
)
# Convertir la columna Barcode al formato barcodeXX en el segundo intento
el_met$barcode <- sprintf("barcode%02d", el_met$barcode)
print(el_met)

#Para MW_P2
# Crear el data frame con los datos proporcionados
mw_met <- data.frame(
  sample = c("AR070_1", "AR070_1", "AR070_1", "AR070_1",
             "AR070_2", "AR070_2", "AR070_2",
             "AR070_3", "AR070_3", "AR070_3", "AR070_3",
             "Mock Community"),
  barcode = c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13)
)
# Añadir las columnas 'HMA' y 'Date' con valores específicos
mw_met$place <- "HMA"
mw_met$date <- "06/02/2024"
mw_met$barcode <- sprintf("barcode%02d", mw_met$barcode)

##AÑADIENDO ULTIMO SECUENCIAMIENTO
last_temp <- data.frame(
  sample = c("AR083_1", "AR083_2", "AR083_3", "AR113_1", "AR113_2", "AR113_3", "AR074", "Mock Community", "Negative Control - Human DNA"),
  barcode = paste0("barcode0",c(1, 2, 3, 4, 5, 6, 7, 8, 9)),
  place = c("HCH", "HCH", "HCH", "INSN", "INSN", "INSN", "HMA", NA, NA),
  date = c("24/04/2024", "24/04/2024", "24/04/2024", "08/07/2024", "08/07/2024", "08/07/2024", "12/03/2024", NA, NA))

seq2111_met<- data.frame(
  sample = c("AR073", "AR078", "AR081", "AR087", "AR088", "AR089", "AR090", "AR092", 
                "Mock Community", "Negative Control - Human DNA"),
  barcode = paste0("barcode",c(73, 74, 75, 76, 77, 78, 79, 80, 81, 82)),
  place = c("HMA", "HMA", "HCH", "HMA", "HCH", "INSN", "HMA", "INSN", NA, NA),
  date = c("05/03/2024", "10/04/2024", "17/04/2024", "08/05/2024", "08/05/2024", "13/05/2024",
           "15/05/2024", "20/05/2024", NA, NA))
seq1212_met <- data.frame(
sample = c("AR097", "AR098", "AR099", "AR100", "AR101", "AR102", "AR103", "AR104", "AR105", "AR106", "Mock Community", "Negative control - Human DNA"),
barcode = paste0("barcode",c(83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94)),
place = c("HCH", "INSN", "HMA", "HCH", "INSN", "HMA", "HCH", "INSN", "HMA", "HCH", NA, NA),
 date = c("29-05-2024", "03-06-2024", "05-06-2024", "05-06-2024", "10-06-2024", "12-06-2024", "12-06-2024", "17-06-2024", "19-06-2024", "19-06-2024", NA, NA)
)
# Unir los `data.frames` de metadata (`el_met`, `temp_met`, `mw_met`) en un único `data.frame`
metadata <- bind_rows(
  "1405Gb" = el_met, 
  final_combined = temp_met, 
  MA_DS9_Crowfoot = mw_met, 
  .id = "source",
  "seq2111"=seq2111_met,
  "DS_3.8GB"=last_temp,
  "seq1212"=seq1212_met
) %>%
  select(source, sample, barcode, place, date) %>% 
  rename(sample_id=sample) %>% 
  mutate(sample = str_c(barcode, source, sep = "_")) %>% 
  select(-c(source))

#En esta sección añadimos el número total de reads para cada barcode, para obtener la abundancia relativa 
last_seq_total<-read_tsv(here("data","countreads_17sset.tsv.tsv"), col_names = c("sample","total"))%>% 
mutate(sample=paste0("17_set_",sample))
temporal_total<-read_tsv(here("data","countreads_ch.tsv.tsv"), col_names = c("sample","total"))
mw_total<-read_tsv(here("data","countreads_mw.tsv.tsv"), col_names = c("sample","total")) 
seq2111_total<-read_tsv(here("data","countreads_seq211.tsv"), col_names = c("sample","total"))%>%
  mutate(sample=paste0(sample,"_seq2111"))
seq1212_total<-read_tsv(here("data","countreads_seq1212.tsv"), col_names = c("sample","total")) %>% mutate(sample=paste0(sample,"_seq1212"))

total_reads<-bind_rows(mw_total,temporal_total,last_seq_total, seq2111_total,seq1212_total)
metadata_dos<-metadata %>% 
  inner_join(total_reads,by="sample")

```
## Procesar resultado de Sylph+ConQuR: EN DESUSO
```{r}
#conq<-read_tsv(here("data","ConQuR_corrected_matrix.tsv"))
# Transponer `conq` usando pivot_longer para crear columnas de 'genome' y 'conteo'
#conq_long <- conq %>% pivot_longer(-sample, names_to = "genome", values_to = "conteo")
#conq_long_sep<- conq_long %>% mutate(sample = str_remove(sample, "\\.fastq"))
# Paso 3: Unir `conq_long` con `metadata` por la columna 'sample'
final_df <- conq_long_sep %>%
  inner_join(metadata_dos, by = "sample") %>% 
  rename(name=genome) %>% 
  rename(sample_conqur=sample) %>% 
  rename(sample=sample_id) %>% 
  rename(directly_reads=conteo) %>% 
  mutate(fraction =directly_reads/total*100) %>% 
  unite("identifier", place, date, sample, sep = "_", remove = FALSE) %>% 
  filter(str_detect(name, bacteria)) %>% 
  select(-c("total","sample_conqur"))

```
## Procesar resultado de sylph
Diferentes inputs
```{r}
#Escoger especies
bacteria=c("Klebsiella pneumoniae", "Pseudomonas aeruginosa",
           "Acinetobacter baumanni", "Escherichia coli", 
           "Enterobacter cloacae","Providencia rettgeri")
#Matriz con conteos de Sylph no corregidos por ConQuR
sylph <- read_tsv(here("data", "abundance_data_export.tsv"))
seq2111 <- read_tsv(here("data","final_data_seq211.tsv"))%>% 
  rename(conteo=count_reads_each_specie,genome=Contig_name,
         sample=Sample_file) %>% 
  mutate(sample=str_replace(sample,"^([^\\.]*)(\\..*)","\\1_seq2111"))
seq1212 <- read_tsv(here("data","final_data_seq1212.tsv"))%>% select(1,4,15) %>% 
  rename(abundance=Sequence_abundance,genome=Contig_name,
         sample=Sample_file) %>% 
  mutate(sample=str_extract(sample,"barcode\\d+"),sample=paste0(sample,"_seq1212")) %>% right_join(seq1212_total, by="sample") %>% mutate(conteo=round((abundance*total)/100)) %>% select(c(1,5,3))

sylph_long<-sylph %>% pivot_longer(
  cols = contains("barcode"), # Selecciona las columnas que contienen los datos de las muestras
  names_to = "sample",          # Nombre de la nueva columna que contendrá los nombres de las muestras
  values_to = "conteo"           # Nombre de la nueva columna que contendrá los valores
) %>% rename(genome=Contig_name)
sylph_long_sep<-sylph_long %>% 
  mutate(sample = str_remove(sample, "\\.fastq")) %>% 
  bind_rows(seq2111,seq1212)
# Paso 3: Unir `conq_long` con `metadata` por la columna 'sample'
final_df <- sylph_long_sep %>%
  inner_join(metadata_dos, by = "sample") %>% 
  rename(name=genome) %>% 
  rename(sample_conqur=sample) %>% 
  rename(sample=sample_id) %>% 
  rename(directly_reads=conteo) %>% 
  filter(map_lgl(name, ~ any(str_detect(.x, bacteria)))) %>% 
  mutate(fraction =directly_reads/total*100,especie = map_chr(name, ~ {
    match <- bacteria[str_detect(.x, bacteria)]
    ifelse(length(match) > 0, match, "Desconocida")
  })) %>% 
  select(-c("sample_conqur"))
print(final_df)
dim(final_df)
print(sum(final_df$directly_reads))
```
## Figura de variación temporal de abundancia de patógenos
```{r}
# Convertir 'fraction' a numérico y 'date' a formato fecha
combined_df2 <- final_df %>%
  mutate(
    fraction = as.numeric(fraction),
    date = dmy(date))%>%  
  filter(
     !str_detect(sample, "Mock ") & 
    !str_detect(sample, "Negative")
   ) %>%
  mutate(sample=str_replace(sample,"(^.*)(_.*)","\\1"))

# Resumir los datos para las líneas (mediana de fracción por fecha, lugar y especie)
data_summary <- combined_df2 %>%
  group_by(date, place, especie) %>%
  summarise(
    median_fraction = median(fraction, na.rm = TRUE),  # Mediana
    q1 = quantile(fraction, 0.25, na.rm = TRUE),  # Cuartil inferior
    q3 = quantile(fraction, 0.75, na.rm = TRUE)   # Cuartil superior
  ) %>%
  ungroup()
gg <- ggplot(combined_df2, aes(x = date, y = fraction, color = place, shape = place)) +
  # Línea de la mediana
  geom_line(
    data = data_summary,
    aes(y = median_fraction, group = place),
    size = 0.8
  ) +
  # Puntos de la mediana con cuartiles (Q1 y Q3)
  geom_pointrange(
    data = data_summary,
    aes(y = median_fraction, ymin = q1, ymax = q3),
    position = position_dodge(width = 0.7)
  ) +
  facet_wrap(~especie, scales = "free_y") +
  labs(
        title = "Percentage abundance of reads assigned to pathogens by date and place",
    x = "Date",
    y = "Relative abundance (%)",
    color = "Place",
    shape = "Place"
  ) +
  theme_minimal() +
  theme(
     plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "italic"),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank()
  )
ggsave("temporal_variation_pathogen.png", plot = gg, width = 8, height = 6, dpi = 300, bg="white")
```

## Grafico de variacion del tamaño de librería
```{r}
combined_librarysize <- metadata_dos %>%
  mutate(
    total = as.numeric(total),
    date = dmy(date))%>%  
  filter(
     !str_detect(sample_id, "Mock ") & 
    !str_detect(sample_id, "Negative")
   ) %>%
  mutate(sample=str_replace(sample_id,"(^.*)(_.*)","\\1")) %>% select(-sample_id)
# Resumir los datos para las líneas (mediana de fracción por fecha, lugar y especie)
data_libsize <- combined_librarysize %>%
  group_by(date, place) %>%
  summarise(
    median_fraction = median(total, na.rm = TRUE),  # Mediana
    q1 = quantile(total, 0.25, na.rm = TRUE),  # Cuartil inferior
    q3 = quantile(total, 0.75, na.rm = TRUE)   # Cuartil superior
  ) %>%
  ungroup()
gg <- ggplot(combined_librarysize, aes(x = date, y = total, color = place, shape = place, fill = place)) +
  # Boxplots por fecha y sample
  geom_boxplot(aes(group = interaction(date, sample)), 
               position = position_dodge(width = 0.7),
               outlier.shape = NA) +
  # Línea de la mediana por lugar
  geom_line(data = data_libsize, 
            aes(x = date, y = median_fraction, group = place),
            size = 1) +
  # Puntos para la mediana
  geom_point(data = data_libsize, 
             aes(x = date, y = median_fraction), 
             size = 1.5) +
  # Etiquetas y tema
  labs(
    title="Total read count per library",
    x = "Date", 
    y = "Count",
    fill = "Place",   # Leyenda para fill
    color = "Place",  # Leyenda para color
    shape = "Place"   # Leyenda para shape
  ) +
  theme_minimal() +
  theme(
     plot.title = element_text(hjust = 0.5,face = "bold"),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank()
  )
ggsave("libsize_temporal_variation_.png", plot = gg, width = 8, height = 6, dpi = 300, bg="white")

```

## Gráfico con adición cruda de data clínica
```{r}

# Segunda tabla: Resultados del análisis de las muestras
date_bacteria <- tibble(
  date = c("04/03/2024", "09/03/2024", "09/03/2024", "12/03/2024", 
                "12/03/2024","12/03/2024", "22/03/2024", "26/03/2024", 
                "30/03/2024", "27/03/2024", "04/09/2024", "07/04/2024", 
                "13/04/2024", "17/04/2024","04/09/2024"),
  especie = c("Escherichia coli", "Klebsiella pneumoniae", "Klebsiella pneumoniae", 
               "Escherichia coli", "Klebsiella pneumoniae", "Klebsiella pneumoniae",
               "Escherichia coli", "Pseudomonas aeruginosa", "Pseudomonas aeruginosa", 
               "Providencia rettgeri", "Pseudomonas aeruginosa", "Klebsiella pneumoniae", 
               "Klebsiella pneumoniae", "Escherichia coli","Pseudomonas aeruginosa")
)

# Convertir la columna DATE_SAMP a formato de fecha
date_bacteria <- date_bacteria %>%
  mutate(date = dmy(date))
gg <- ggplot(combined_df2, aes(x = date, y = fraction, color = place, shape = place)) +
  # Línea de la mediana
  geom_line(
    data = data_summary,
    aes(y = median_fraction, group = place),
    size = 0.8
  ) +
  # Puntos de la mediana con cuartiles (Q1 y Q3)
  geom_pointrange(
    data = data_summary,
    aes(y = median_fraction, ymin = q1, ymax = q3),
    position = position_dodge(width = 0.7)
  ) +
  facet_wrap(~especie, scales = "free_y") +
  labs(
        title = "Percentage abundance of reads assigned to pathogens by date and place",
    x = "Date",
    y = "Relative abundance (%)",
    color = "Place",
    shape = "Place"
  ) +
    # Aquí colocamos los puntos en y=0 o algún valor bajo, sólo para indicar la presencia.
  geom_point(
    data = date_bacteria,
    aes(x = date, y = 0),
    color = "black",
    shape = 16,
    fill = "black",
    size = 2,
    inherit.aes = FALSE,
  ) +
  theme_minimal() +
  theme(
     plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "italic"),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank()
  )
ggsave("temporal_variation_pathogen_clinical.png", plot = gg, width = 8, height = 6, dpi = 300)

```
## Correlaciones
```{r}
# Gráfico 1: Abundance vs Total Reads
p1 <- ggplot(final_df, aes(x = total, y = fraction)) +
  geom_point(alpha = 0.6, color = "darkorange") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red") +
  stat_cor(
    aes(label = paste0(after_stat(rr.label), "~", after_stat(p.label))),
    color = "black",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 3,
    parse = TRUE
  ) +
  facet_wrap(~ especie, scales = "free") +
  labs(
    title = "Correlation between relative abundance and total reads",
    x = "Total reads",
    y = "Relative abundance (%)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "italic"),  # Titles of facets (italicized species)
    plot.title = element_text(hjust = 0.5), # Centered italic title
    axis.text.x = element_text(size = 10), # Reduced axis text size
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11) # Reduced axis title size
  )

# Gráfico 2: Direct Reads vs Total Reads
p2 <- ggplot(final_df, aes(x = total, y = directly_reads)) +
  geom_point(alpha = 0.6, color = "darkorange") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "red") +
  stat_cor(
    aes(label = paste0(after_stat(rr.label), "~", after_stat(p.label))),
    color = "black",
    label.x.npc = "left",
    label.y.npc = "top",
    size = 3,
    parse = TRUE
  ) +
  facet_wrap(~ especie, scales = "free") +
  labs(
    title = "Correlation between assigned reads and total reads",
    x = "Total reads",
    y = "Assigned reads"
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text = element_text(face = "italic"),  # Titles of facets (italicized species)
    plot.title = element_text(hjust = 0.5), # Centered italic title
    axis.text.x = element_text(size = 10), # Reduced axis text size
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 11) # Reduced axis title size
  )

# Mostrar gráficos
ggsave("correlation_abundancevstotalreads.png", plot = p1, width = 8, height = 6, dpi = 300, bg="white")
ggsave("correlation_assignedreadsvstotalreads.png", plot = p2, width = 8, height = 6, dpi = 300, bg="white")
```

