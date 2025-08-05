---
title: "exploratoyr_pathogen_analysis"
output: html_document
date: "2024-10-04"
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## What is this document about?

Performs exploratory analysis for abundance metagenomic data variation accross time and sampling sites.
Per-sample bars, single species, method comparison (Sylph/ConQuR/Kraken), Mock check. Please take a look at the produced images in the `exploratory_pathogen_analysis` folder

## Let's load libraries and define color schemes
```{r}
library(tidyverse)
library(readr)
library(here)
library(janitor)
library(ggplot2)
library(scales)
library(boot)
library(ggh4x)
library(RColorBrewer)
setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/loui/abundance")
```

```{r}
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

## Upload abundance data from Kraken2+Bracken
```{r}
ael<-read_tsv("all_el.tsv.tsv")
atemp<-read_tsv("all_temp.tsv.tsv")
amw<-read_tsv("combined_summary.tsv.txt")
colnames(ael) <- tolower(colnames(ael))
colnames(amw) <- tolower(colnames(amw))
colnames(atemp) <- tolower(colnames(atemp))
```

## Prepare the metadata for each barcode from each sequencing folder.
```{r}
# Create the data frame with the provided data
met_temp <- data.frame(
  `sample` = c("AR071_1", "AR071_2", "AR071_2", "AR071_2", "AR071_3",
     "AR072_1", "AR072_1", "AR072_1", "AR072_2", "AR072_2",
     "AR072_2", "AR072_3", "Mock Community", "Negative Control - Human DNA"),
  barcode = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
)

# Convert the Barcode column to the barcodeXX format
met_temp$barcode <- sprintf("barcode%02d", met_temp$barcode)

# Show the resulting data frame
print(met_temp)

# Create the second data frame with the new data
additional_data <- data.frame(
  sample = c("AR071_1", "AR071_2", "AR071_3", 
     "AR072_1", "AR072_2", "AR072_3"),
  place = rep("HMA", 6),  # 'HMA' column with constant value
  date = c("13/02/2024", "13/02/2024", "13/02/2024", 
     "20/02/2024", "20/02/2024", "20/02/2024")
)
# Merge both data frames by the 'sample_of_extraction' column
temp_met <- merge(met_temp, additional_data, by = "sample", all.x = TRUE)

# For L folder
# Create the third data frame with the data from the 2nd attempt
el_met<- data.frame(
  sample = c("AR086", "AR093","AR093","AR093", "AR095", 
     "Mock Community", "Negative Control - Human DNA"),
  `date` = c("06/05/2024", "22/05/2024", "22/05/2024", "22/05/2024", "27/05/2024", NA, NA),
  place = c("INSN", "HMA", "HMA", "HMA", "INSN",NA, NA),
  `barcode` = c(24, 23, 22, 21, 20, 19, 18)
)
# Convert the Barcode column to the barcodeXX format in the second attempt
el_met$barcode <- sprintf("barcode%02d", el_met$barcode)
print(el_met)

# For MW_P2
# Create the data frame with the provided data
mw_met <- data.frame(
  sample = c("AR070_1", "AR070_1", "AR070_1", "AR070_1",
     "AR070_2", "AR070_2", "AR070_2",
     "AR070_3", "AR070_3", "AR070_3", "AR070_3",
     "Mock Community"),
  barcode = c(1, 2, 3, 4, 5, 6, 7, 9, 10, 11, 12, 13)
)
# Add the 'HMA' and 'Date' columns with specific values
mw_met$place <- "HMA"
mw_met$date <- "06/02/2024"
mw_met$barcode <- sprintf("barcode%02d", mw_met$barcode)

##Adding last sequencing metadata
last_temp <- data.frame(
  sample_id = c("AR083_1", "AR083_2", "AR083_3", "AR113_1", "AR113_2", "AR113_3", "AR074", "Mock Community", "Negative Control - Human DNA"),
  barcode = paste0("barcode0",c(1, 2, 3, 4, 5, 6, 7, 8, 9)),
  place = c("HCH", "HCH", "HCH", "INSN", "INSN", "INSN", "HMA", NA, NA),
  date = c("24/04/2024", "24/04/2024", "24/04/2024", "08/07/2024", "08/07/2024", "08/07/2024", "12/03/2024", NA, NA))

seq2111_met<- data.frame(
  sample_id = c("AR073", "AR078", "AR081", "AR087", "AR088", "AR089", "AR090", "AR092", 
    "Mock community", "Negative Control - Human DNA"),
  barcode = paste0("barcode",c(73, 74, 75, 76, 77, 78, 79, 80, 81, 82)),
  place = c("HMA", "HMA", "HCH", "HMA", "HCH", "INSN", "HMA", "INSN", NA, NA),
  date = c("05/03/2024", "10/04/2024", "17/04/2024", "08/05/2024", "08/05/2024", "13/05/2024",
     "15/05/2024", "20/05/2024", NA, NA))

# Merge the metadata data.frames (`el_met`, `temp_met`, `mw_met`) into a single `data.frame`
metadata <- bind_rows(
  "1405Gb" = el_met, 
  final_combined = temp_met, 
  MA_DS9_Crowfoot = mw_met, 
  .id = "source"
) %>%
  select(source, sample, barcode, place, date) %>% 
  rename(sample_id=sample) %>% 
  mutate(sample = str_c(barcode, source, sep = "_")) %>% 
  select(-c(source))
last_temp<-last_temp %>% mutate(sample=str_c("17_set",barcode, sep="_"))

metadata<-bind_rows(metadata,last_temp, seq2111_met%>%mutate(sample=paste0(barcode,".seq2111")))
print(metadata)
```

## Keep metadata processing and join metadata with abundance data from each main folder. At the end of this chunck, we filter to get the info from the abundance of one pathogen of interest
```{r}
# Transpose and restructure the data frame
atemp_long <- atemp %>%
  # Select the species columns and the barcode columns that contain reads and fraction
  select(name, starts_with("barcode")) %>%
  # Convert the reads and fraction columns to long format
  pivot_longer(cols = -name, 
     names_to = c("barcode", "type"), 
     names_pattern = "(barcode\\d+).(.+)") %>%
  # Rearrange so that each row is a barcode-species with its reads and fraction values
  pivot_wider(names_from = type, values_from = value)

# Transpose and restructure the data frame
ael_long <- ael %>%
  # Select the species columns and the barcode columns that contain reads and fraction
  select(name, starts_with("barcode")) %>%
  # Convert the reads and fraction columns to long format
  pivot_longer(cols = -name, 
     names_to = c("barcode", "type"), 
     names_pattern = "(barcode\\d+).(.+)") %>%
  # Rearrange so that each row is a barcode-species with its reads and fraction values
  pivot_wider(names_from = type, values_from = value)
# Show the final table with columns: "barcode", "name", "reads", "fraction"
atemp_long
ael_long

mw_of<-merge(amw, mw_met, by = "barcode", all.x = TRUE)
temp_of<-merge(atemp_long, temp_met, by = "barcode", all.x = TRUE)
el_of<-merge(ael_long, el_met, by = "barcode", all.x = TRUE)

mw_of <- mw_of %>%
  select(-kraken_assigned_reads, -added_reads,-taxonomy_id,-taxonomy_level)%>% 
  rename(directly_reads=new_estimated_reads) %>%
  rename(fraction=fraction_total_reads) %>% 
  mutate(fraction=fraction*100)

# 2. Remove the 'total_reads' column from the three data frames
el_of <- el_of %>%
  rename(fraction=cladereads..) %>% 
  rename(directly_reads=cladereads)

temp_of <- temp_of %>%
  select(-starts_with("taxonre")) %>% 
  rename(fraction=cladereads..) %>% 
  rename(directly_reads=cladereads)

el_of<- el_of%>%
  mutate(fraction = as.numeric(fraction))
temp_of<- temp_of %>% 
  mutate(fraction = as.numeric(fraction))
el_of2 <- el_of %>% filter(!is.na(date))
temp_of2 <- temp_of %>% filter(!is.na(date))
mw_of2 <- mw_of %>% filter(!is.na(date)) %>% filter(sample!="Mock Community")
##DEFIne species
esp<-"Escherichia coli"
###Filter by species
el_ofs<-el_of2 %>% filter(name==esp)
temps_ofs<-temp_of2 %>% filter(name==esp)
mw_ofs<-mw_of2 %>% filter(name==esp)
```

## Calculate fraction and read count for the pathogen in each sample (barcode) and make graph
```{r}
# Join the columns 'place', 'date', and 'sample' into 'identifier' without removing the original columns
el_off <- el_ofs %>% unite("identifier", place, date, sample, sep = "_", remove = FALSE)
temp_off <- temps_ofs %>% unite("identifier", place, date, sample, sep = "_", remove = FALSE)
mw_off <- mw_ofs %>% unite("identifier", place, date, sample, sep = "_", remove = FALSE)

# Join the three data frames
combined_df <- bind_rows(el_off, temp_off, mw_off)

# Convert 'fraction' to numeric and 'date' to date format
combined_df <- combined_df %>%
  mutate(
  fraction = as.numeric(fraction),
  date = dmy(date)
  )

# Modify 'identifier' to remove 'AR0...' but keep the suffix
combined_df <- combined_df %>%
  mutate(identifier = gsub("AR0\\d+(_\\d+)?", "\\1", identifier)) %>%
  mutate(identifier = gsub("_+", "_", identifier)) %>%
  mutate(identifier = gsub("^_|_$", "", identifier)) %>%
  mutate(identifier = gsub("^[[:alnum:]]+_", "", identifier))

# Create the 'group' column based on 'place'
combined_df <- combined_df %>%
  mutate(group = ifelse(place == "INSN", "INSN", "HMA"))

# Sort 'identifier' by 'group' and 'date'
combined_df <- combined_df %>%
  arrange(factor(group, levels = c("INSN", "HMA")), date)

# Update the factor levels of 'identifier' in the desired order
combined_df$identifier <- factor(combined_df$identifier, levels = unique(combined_df$identifier))

# Calculate 'mean_fraction' and 'mean_directly_reads' by 'identifier'
mean_fraction_df <- combined_df %>%
  group_by(identifier, group) %>%
  summarise(
  mean_fraction = mean(fraction),
  mean_directly_reads = mean(directly_reads),
  .groups = 'drop'
  )

# Ensure that the factor levels of 'identifier' match
mean_fraction_df$identifier <- factor(mean_fraction_df$identifier, levels = levels(combined_df$identifier))

# Create the plot
p <- ggplot() +
  geom_bar(
  data = mean_fraction_df,
  aes(x = identifier, y = mean_fraction, fill = mean_fraction),
  stat = "identity",
  width = 0.6
  ) +
  geom_text(
  data = mean_fraction_df,
  aes(
    x = identifier,
    y = mean_fraction,
    label = paste0("     n=", round(mean_directly_reads))
  ),
  hjust = 0,
  size = 2.5
  ) +
  geom_point(
  data = combined_df,
  aes(x = identifier, y = fraction),
  color = "black",
  size = 1.5,
  position = position_jitter(width = 0.1, height = 0)
  ) +
  scale_fill_gradient(
  low = "#FF7F0E",
  high = "#1F77B4",
  name = "Mean\nPercentage"
  ) +
  coord_flip() +
  # Use facet_wrap2 with free scales on the Y axis
  facet_wrap2(
  ~ group,
  scales = "free",
  ncol = 1
  ) +
  # Adjust the Y axis for each facet
  facetted_pos_scales(
  y = list(
    group == "INSN" ~ scale_y_continuous(labels = function(breaks) {
  # Customize labels: 0.25 with 2 decimals, the rest with 1 decimal
  sapply(breaks, function(x) sprintf("%.1f", x))},
  expand = expansion(mult = c(0, 0.10))),
      group == "HMA" ~ scale_y_continuous(  breaks = if (round(max(combined_df$fraction)) > 3) {
  # If the maximum rounded value is greater than 3, only show integers
  seq(1, round(max(combined_df$fraction)), by = 1)
} else {
  # Otherwise, show the specified values
  c(0, 0.25, 0.5, seq(1, round(max(combined_df$fraction))))
},labels = function(breaks) {
  # Customize labels: 0.25 with 2 decimals, the rest with 1 decimal
  sapply(breaks, function(x) ifelse(x == 0.25, sprintf("%.2f", x), sprintf("%.1f", x)))
  },
expand = expansion(mult = c(0, 0.10)))))
  +
  theme_minimal(base_size = 14) +
  labs(
  title = bquote("Mean percentage of " * italic(.(esp)) * " reads"),
  x = "Sample",
  y = "Read percentage"
  ) +
  theme(
  axis.text.y = element_text(size = 12),
  axis.title.y = element_text(size = 14),
  axis.title.x = element_text(size = 14),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.position = "right",
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
  panel.grid.minor = element_blank()
  )

# Show the plot
print(p)
```

## Process Sylph+CONQUR abundance normalization results. 
```{r}
conq<-read_tsv("ConQuR_corrected_matrix.tsv")
sylph<-read_tsv("abundance_data_export.tsv")
library(tidyverse)

# Step 1: Transpose `conq` using pivot_longer to create 'genome' and 'conteo' columns
conq_long <- conq %>% pivot_longer(-sample, names_to = "genome", values_to = "conteo")
sylph_long<-sylph %>% pivot_longer(
    cols = contains("barcode"), # Selects the columns containing the sample data
  names_to = "sample",          # Name of the new column that will contain the sample names
  values_to = "conteo"           # Name of the new column that will contain the values
  ) %>% rename(genome=Contig_name)
##ADDING LAST SEQUENCING
last_temp <- data.frame(
  sample_id = c("AR083_1", "AR083_2", "AR083_3", "AR113_1", "AR113_2", "AR113_3", "AR074", "Mock Community", "Negative Control - Human DNA"),
  barcode = paste0("barcode0",c(1, 2, 3, 4, 5, 6, 7, 8, 9)),
  place = c("HCH", "HCH", "HCH", "INSN", "INSN", "INSN", "HMA", NA, NA),
  date = c("24/04/2024", "24/04/2024", "24/04/2024", "08/07/2024", "08/07/2024", "08/07/2024", "12/03/2024", NA, NA))
# Step 2: Join the metadata data.frames (`el_met`, `temp_met`, `mw_met`) into a single `data.frame`


seq2111_met<- data.frame(
  sample_id = c("AR073", "AR078", "AR081", "AR087", "AR088", "AR089", "AR090", "AR092", 
        "Mock community", "Negative Control - Human DNA"),
  barcode = paste0("barcode",c(73, 74, 75, 76, 77, 78, 79, 80, 81, 82)),
  place = c("HMA", "HMA", "HCH", "HMA", "HCH", "INSN", "HMA", "INSN", NA, NA),
  date = c("05/03/2024", "10/04/2024", "17/04/2024", "08/05/2024", "08/05/2024", "13/05/2024",
       "15/05/2024", "20/05/2024", NA, NA))

metadata<-bind_rows(metadata,last_temp, seq2111_met%>%mutate(sample=barcode))


metadata <- bind_rows(
  "1405Gb" = el_met, 
  final_combined = temp_met, 
  MA_DS9_Crowfoot = mw_met, 
  .id = "source"
) %>%
  select(source, sample, barcode, place, date) %>% 
  rename(sample_id=sample) %>% 
  mutate(sample = str_c(barcode, source, sep = "_")) %>% 
  select(-c(source))
last_temp<-last_temp %>% mutate(sample=str_c("17_set",barcode, sep="_"))
metadata<-bind_rows(metadata,last_temp)

##HERE ADD THE TOTAL NUMBER OF READS TO FIND THE FRACTION OF BAUMANNII
last_seq_total<-read_tsv("countreads_17sset.tsv.tsv", col_names = c("sample","total"))%>% 
  mutate(sample=paste0("17_set_",sample))
temporal_total<-read_tsv("countreads_ch.tsv.tsv", col_names = c("sample","total"))
mw_total<-read_tsv("countreads_mw.tsv.tsv", col_names = c("sample","total")) 
seq2111_total<-read_tsv("countreads_seq211.tsv", col_names = c("sample","total"))

total_reads<-bind_rows(mw_total,temporal_total,last_seq_total, seq2111_total)

sylph_long_sep<-sylph_long %>% 
  mutate(sample = str_remove(sample, "\\.fastq"))
conq_long_sep<- conq_long %>%
  mutate(sample = str_remove(sample, "\\.fastq"))
metadata_dos<-metadata %>% 
  inner_join(total_reads,by="sample")

bacteria="Klebsiella pneumoniae"
# Step 3: Join `conq_long` with `metadata` by the 'sample' column
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
## Now let's see how the figure changes for the same pathogen using these tools

```{r}
# Convert 'fraction' to numeric and 'date' to date format
combined_df2 <- final_df %>%
  mutate(
  fraction = as.numeric(fraction),
  date = dmy(date))# %>%  
  #filter(
   #!str_detect(sample, "Mock Community") & 
  # !str_detect(sample, "Negative Control")
   #)
combined_df2 <- combined_df2 %>%
  separate(identifier, into = c("parte1", "fecha", "parte3", "sufijo"), sep = "_") %>%
  unite("identifier", fecha, sufijo, sep = "_") %>% select(-c("parte1","parte3")) %>% mutate(identifier = str_remove(identifier, "_NA$"))

# Make sure the data is ordered as before

# Create the 'group' column based on 'place' (as before)
combined_df2 <- combined_df2 %>%
  mutate(group = case_when(
  sample =="Mock Community" ~"Mock Community",
  str_detect(sample, "Negative Control") ~ "Negative Control",
  place == "INSN" ~ "INSN",
  place == "HMA" ~ "HMA",
  TRUE ~ "HCH"  # Assigns "HCH" if it is neither "INSN" nor "HMA"
  ))

# Extract the date and number from 'identifier'
combined_df2 <- combined_df2 %>%
  mutate(
  number_from_id = as.numeric(sub(".*_", "", identifier))
  )

# Modify 'identifier' only for HMA
combined_df2 <- combined_df2 %>%
  mutate(
  identifier = case_when(
    sample == "Mock Community" ~ "Mock Community",
    str_detect(sample, "Negative Control") ~ "Negative Control",
    place == "HMA" & date == as.Date("2024-02-06") ~ "HMA_baseline",
    place == "HMA" & date == as.Date("2024-02-13") ~ "HMA_1_week",
    place == "HMA" & date == as.Date("2024-02-20") ~ "HMA_2_weeks",
    place == "HMA" & date == as.Date("2024-03-12") ~ "HMA_5_weeks",
    place == "HMA" & date == as.Date("2024-05-22") ~ "HMA_3_months",
    place == "INSN" & date == as.Date("2024-05-06") ~ "INSN-I_I",
    place == "INSN" & date == as.Date("2024-05-27")~
    "INSN-I_II",
    place == "INSN" & date == as.Date("2024-07-08") ~ "INSN-II",
    TRUE ~ place
  ))

palette_colors <- c(
  "HMA_baseline"   = "#2171B5",  # Darker blue
  "HMA_1_week"     = "#4292C6",  # Intermediate blue
  "HMA_2_weeks"    = "#6BAED6",  # Lighter blue
  "HMA_5_weeks"    = "#9ECAE1",  # Even lighter blue
  "HMA_3_months"   = "#2131B5",  # Custom blue
  "INSN-I_I"         = "#74C476",  # Intermediate green
  "INSN-I_II"      = "#74C476",
  "INSN-II"        = "#A1D99B",  # Lighter green
  "HCH"            = "#6A51A3",  # Dark purple
  "Mock Community" = "#F16913",  # Medium orange
  "Negative Control" = "#CB181D" # Dark red
)
# Sort 'identifier' by 'group', 'date', and 'number_from_id'
combined_df2 <- combined_df2 %>%
  arrange(
  factor(group, levels = c("INSN", "HCH", "HMA")),
  date,
  number_from_id
  )

# Update the factor levels of 'identifier' in the desired order
combined_df2$identifier <- factor(combined_df2$identifier, levels = unique(combined_df2$identifier))

# Calculate 'mean_fraction' and 'mean_directly_reads' by 'identifier' and 'group'
mean_fraction_df2 <- combined_df2 %>%
  group_by(identifier, group) %>%
  summarise(
  mean_fraction = mean(fraction),
  mean_directly_reads = mean(directly_reads),
  .groups = 'drop'
  )

# Ensure that the factor levels of 'identifier' match in both dataframes
mean_fraction_df2$identifier <- factor(mean_fraction_df2$identifier, levels = levels(combined_df2$identifier))

# Split the data into two dataframes: one for HMA and another for INSN and HCH
combined_df2_hma <- combined_df2 %>% filter(group %in% c("HMA", "Mock Community", "Negative Control"))
mean_fraction_df2_hma <- mean_fraction_df2 %>% filter(group %in% c("HMA", "Mock Community", "Negative Control"))

combined_df2_others <- combined_df2 %>% filter(group %in% c("HCH", "INSN", "Mock Community", "Negative Control"))
mean_fraction_df2_others <- mean_fraction_df2 %>% filter(group %in% c("HCH", "INSN", "Mock Community", "Negative Control"))
# Create the plot for HMA and store it in 'p_hma'
p_hma <- ggplot() +
  geom_bar(
  data = mean_fraction_df2_hma,
  aes(x = identifier, y = mean_fraction, fill = identifier),
  stat = "identity",
  width = 0.6,
  show.legend=FALSE
  ) +
  geom_text(
  data = mean_fraction_df2_hma,
  aes(
    x = identifier,
    y = mean_fraction,
    label = paste0("     n=", round(mean_directly_reads))
  ),
  hjust = 0,
  size = 2.5
  ) +
  geom_point(
  data = combined_df2_hma,
  aes(x = identifier, y = fraction),
  color = "black",
  size = 1.5,
  position = position_jitter(width = 0.1, height = 0)
  ) +
  scale_fill_manual(
  values = palette_colors,
  name = "Identifier"
  ) +
  coord_flip() +
  # Adjust the Y axis for HMA
  scale_y_continuous(
  breaks = if (round(max(combined_df2_hma$fraction)) > 3) {
    seq(1, round(max(combined_df2_hma$fraction)), by = 1)
  } else {
    c(0, 0.25, 0.5, seq(1, round(max(combined_df2_hma$fraction))))
  },
  labels = function(breaks) {
    sapply(breaks, function(x) ifelse(x == 0.25, sprintf("%.2f", x), sprintf("%.1f", x)))
  },
  expand = expansion(mult = c(0, 0.10))
  ) +
  theme_minimal(base_size = 14) +
  labs(
  title = bquote("Mean percentage of " * italic(.(bacteria)) * " reads - HMA"),
  x = "Sample",
  y = "Read percentage"
  ) +
  theme(
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 14),
  axis.title.x = element_text(size = 14),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.position = "right",
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
  panel.grid.minor = element_blank()
  )

# Create the plot for INSN and HCH and store it in 'p_others'
p_others <- ggplot() +
  geom_bar(
  data = mean_fraction_df2_others,
  aes(x = identifier, y = mean_fraction, fill = identifier),
  stat = "identity",
  width = 0.6,
  show.legend=FALSE
  ) +
  geom_text(
  data = mean_fraction_df2_others,
  aes(
    x = identifier,
    y = mean_fraction,
    label = paste0("     n=", round(mean_directly_reads))
  ),
  hjust = 0,
  size = 2.5
  ) +
  geom_point(
  data = combined_df2_others,
  aes(x = identifier, y = fraction),
  color = "black",
  size = 1.5,
  position = position_jitter(width = 0.1, height = 0)
  ) +
  scale_fill_manual(
  values = palette_colors,
  name = "Identifier"
  ) +
  coord_flip() +
  # Adjust the Y axis for INSN and HCH
  scale_y_continuous(
  labels = function(breaks) {
    sapply(breaks, function(x) sprintf("%.1f", x))
  },
  expand = expansion(mult = c(0, 0.10))
  ) +
  theme_minimal(base_size = 14) +
  labs(
  title = bquote("Mean percentage of " * italic(.(bacteria)) * " reads - INSN and HCH"),
  x = "Sample",
  y = "Read percentage"
  ) +
  theme(
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 14),
  axis.title.x = element_text(size = 14),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 12),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.position = "right",
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
  panel.grid.minor = element_blank()
  )

# Show the plots
print(p_hma)
print(p_others)

```

## Now compare the taxonomic profiling of CONQUR+Sylph with that of Kraken+Bracken in the Mock Community
```{r}
combined_dfk <- bind_rows(el_of, temp_of, mw_of)
combined_dfk<- combined_dfk %>% mutate(name = str_replace(name, "Bacillus subtilis", "Bacillus_not subtilis"))%>% mutate(name = str_replace(name, "Bacillus spizizenii", "Bacillus subtilis"))

bacterias=c(
  "Listeria monocytogenes",
  "Pseudomonas aeruginosa",
  "Bacillus subtilis",
  "Escherichia coli",
  "Salmonella enterica",
  "Lactobacillus fermentum",
  "Enterococcus faecalis",
  "Staphylococcus aureus",
  "Saccharomyces cerevisiae",
  "Cryptococcus neoformans"
)
case_when_expression <- paste0(
  "case_when(",
  paste(
  sprintf("str_detect(name, '%s') ~ '%s'", bacterias, bacterias),
  collapse = ", "
  ),
  ", TRUE ~ NA_character_)"
)
# Step 3: Join `conq_long` with `metadata` by the 'sample' column

two<-conq_long_sep %>%
  mutate(source = "Sylph+ConQuR") %>% 
  bind_rows(sylph_long_sep %>% mutate(source="Sylph"))


final_df <- two %>%
  inner_join(metadata_dos, by = "sample") %>% 
  rename(name=genome) %>% 
  rename(sample_conqur=sample) %>% 
  rename(sample=sample_id) %>% 
  rename(directly_reads=conteo) %>% 
  mutate(fraction =directly_reads/total*100) %>% 
  select(-c("total","sample_conqur")) %>% 
  bind_rows(combined_dfk %>% mutate(source = "Kraken+Bracken")) %>% 
  filter(str_detect(name, paste(bacterias, collapse = "|"))) %>% 
  mutate(species = eval(parse(text = case_when_expression)))

# Convert 'fraction' to numeric and 'date' to date format
combined_df2 <- final_df %>%
  mutate(
  fraction = as.numeric(fraction),
  date = dmy(date))
# Create the 'group' column based on 'place' (as before)
combined_df2 <- combined_df2 %>%
  mutate(group = case_when(
  sample =="Mock Community" ~"Mock Community",
  str_detect(sample, "Negative Control") ~ "Negative Control",
  place == "INSN" ~ "INSN",
  place == "HMA" ~ "HMA",
  TRUE ~ "HCH"  # Assigns "HCH" if it is neither "INSN" nor "HMA"
  )) %>% filter(group %in% c("Mock Community")) %>% 
  filter(!str_detect(name, "Escherichia coli 1.2741 Contig1651129016"))

# Sort 'identifier' by 'group', 'date', and 'number_from_id'
combined_df2 <- combined_df2 %>%
  arrange(
  source,
  factor(group, levels = c("Mock Community", "Negative Control")),
  
  ) %>% 
  mutate(species = str_replace_all(species, " ", "\n"))

# Calculate 'mean_fraction' and 'mean_directly_reads' by 'identifier' and 'group'
mean_fraction_df2 <- combined_df2 %>%
  group_by(species, group, source) %>%
  summarise(
  mean_fraction = median(fraction),
  mean_directly_reads = median(directly_reads),
  .groups = 'drop'
  ) %>%
  arrange(desc(mean_fraction)) %>%
  mutate(species = factor(species, levels = unique(species)))

# Order 'species' in 'combined_df2' based on the order of 'mean_fraction_df2'
combined_df2 <- combined_df2 %>%
  mutate(species = factor(species, levels = levels(mean_fraction_df2$species)))


palette_colors <- brewer.pal(n = 9, name = "Set1")  # You can adjust `n` to the number of groups you have

# Create the improved plot for publication
p_mock <- ggplot() +
  geom_bar(
  data = mean_fraction_df2,
  aes(x = species, y = mean_fraction, fill = species),
  stat = "identity",
  width = 0.6,
  show.legend = FALSE
  ) +
  geom_text(
  data = mean_fraction_df2,
  aes(
    x = species,
    y = mean_fraction,
    label = paste0("n = ", round(mean_directly_reads))
  ),
  hjust=0.5,  # Move slightly to the left
  vjust=-0.5,
  size = 2.8,  # Larger text size for visibility
  color = "black"
  ) +
  geom_point(
  data = combined_df2,
  aes(x = species, y = fraction),
  color = "black",
  size = 1.8,
  shape=16,
  alpha=0.6,# Increase the size of the points for better visibility
  position = position_jitter(width = 0.1, height = 0)
  ) +
  geom_hline(yintercept = 12, color = "red", size = 0.8, alpha=0.7)+
  facet_wrap(~source) +  # Facets in one column for clarity   
  scale_fill_manual(
  values = palette_colors,
  name = "Group"
  )+
  theme_minimal(base_size = 16) +  # Larger base size
  labs(
  title = "Percentage of Reads from Mock Community Species",
  x = "Species",
  y = "Read Percentage"
  ) +
  theme(
  axis.text.y = element_text(size = 10),
  axis.title.y = element_text(size = 12, face = "bold"),
  axis.title.x = element_text(size = 12, face = "bold"),
  axis.text.x = element_text(size = 10, face = "italic", angle = 45, hjust = 1),
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.title = element_text(size = 14, face = "bold"),
  legend.text = element_text(size = 12),
  panel.grid.major = element_line(size = 0.5, linetype = 'solid', colour = "gray90"),
  panel.grid.minor = element_blank(),
  legend.position = "right"
  ) +
  coord_cartesian(ylim = c(0, 18)) 

p_mock

```

