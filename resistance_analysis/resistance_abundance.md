---
title: "Resistance Abundance Analysis"
output: html_document
date: "2024-11-10"
---

## What is this notebook about?

In this notebook I merged per-sample resistance-gene counts (MW/MA_DS9_Crowfoot, L/1405Gb, wastewater/5.2Gb, LIG_P2/3.8Gb), and also with metadata (barcode→sample_id, place, date, batch). 
Then parsed features into Gene and Drug_Class, and joining total read counts to compute percent abundance. 
From that this notebook first filters to carbapenem genes to produce barplots (means per timepoint/identifier with jittered replicates) for HMA and for INSN/HCH, then expands to six Drug_Class categories (penam, tetracycline, fluoroquinolone, cephalosporin, aminoglycoside, carbapenem) to draw temporal trends by place using median lines with IQR ribbons. 
Finally, it overlays clinical data by converting Excel dates, filtering isolates with any R call, and adding a histogram on a secondary axis to compare clinical resistance events against metagenomic trends. 

```{r}

library(tidyverse)

library(readr)
library(here)
library(janitor)
library(ggplot2)
library(scales)
library(boot)
library(ggh4x)
library(openxlsx)

mw_data<-read_tsv(here("data","count_mw_allDS.tsv"))
l_data<-read_tsv(here("data","count_L_1405DS.tsv"))
t_data<-read_tsv(here("data","count_wastewatersch_5.2DS.tsv"))
last<-read_tsv(here("data","count_LIG_P2_3.8DS.tsv"))
```
## PREPARE METADATA AND LIBRARY SIZE

```{r}
met_temp <- data.frame(
  `sample` = c("AR071_1", "AR071_2", "AR071_2", "AR071_2", "AR071_3",
                             "AR072_1", "AR072_1", "AR072_1", "AR072_2", "AR072_2",
                             "AR072_2", "AR072_3", "Mock Community", "Negative control - Human DNA"),
  barcode = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)
)
# Convert the Barcode column to barcodeXX format
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


# Create the third data frame with the 2nd attempt data
 el_met<- data.frame(
  sample = c("AR086", "AR093","AR093","AR093", "AR095", 
                             "Mock Community", "Negative control - Human DNA"),
  `date` = c("06/05/2024", "22/05/2024", "22/05/2024", "22/05/2024", "27/05/2024", NA, NA),
  place = c("INSN", "HMA", "HMA", "HMA", "INSN",NA, NA),
  `barcode` = c(24, 23, 22, 21, 20, 19, 18)
)

# Convert the Barcode column to barcodeXX format in the second attempt
el_met$barcode <- sprintf("barcode%02d", el_met$barcode)

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


last_temp <- data.frame(
  sample_id = c("AR083_1", "AR083_2", "AR083_3", "AR113_1", "AR113_2", "AR113_3", "AR074", "Mock Community", "Negative control - Human DNA"),
  barcode = paste0("barcode0",c(1, 2, 3, 4, 5, 6, 7, 8, 9)),
  place = c("HCH", "HCH", "HCH", "INSN", "INSN", "INSN", "HMA", NA, NA),
  date = c("24/04/2024", "24/04/2024", "24/04/2024", "08/07/2024", "08/07/2024", "08/07/2024", "12/03/2024", NA, NA))

seq1212_met <- data.frame(
  sample_id = c("AR097", "AR098", "AR099", "AR100", "AR101", "AR102", "AR103", "AR104", "AR105", "AR106", "Mock Community", "Negative control - Human DNA"),
  barcode = paste0("barcode",c(83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94)),
  place = c("HCH", "INSN", "HMA", "HCH", "INSN", "HMA", "HCH", "INSN", "HMA", "HCH", NA, NA),
  date = c("29-05-2024", "03-06-2024", "05-06-2024", "05-06-2024", "10-06-2024", "12-06-2024", "12-06-2024", "17-06-2024", "19-06-2024", "19-06-2024", NA, NA)
)

el_met<-el_met %>% mutate(folder="L_folder")
temp_met<-temp_met %>% mutate(folder="temporal")
mw_met<-mw_met %>% mutate(folder="MW_P2")
last_temp<-last_temp %>% mutate(folder="LIG_P2")

metadata <- bind_rows(
  "1405Gb" = el_met, 
  final_combined = temp_met, 
  MA_DS9_Crowfoot = mw_met, 
  .id = "source"
) %>%
  select(source, sample, barcode, place, date, folder) %>% 
  rename(sample_id=sample) %>% 
  mutate(sample = str_c(barcode, source, sep = "_")) %>% 
  select(-c(source))

last_temp<-last_temp %>% mutate(sample = paste0(barcode,"_DS_3.8GB"))
seq1212_met<-seq1212_met %>% mutate(sample=paste0(barcode,"_seq1212"))
metadata<-bind_rows(metadata,last_temp, seq1212_met) %>% 
  filter(!str_detect(sample_id, "Negative"))

last_seqDS<-read.table(here("data","countreadsDS_17sset.tsv.tsv"), col.names = c("sample","total"), sep= " ")%>% mutate(sample = gsub(".fastq.gz.*$","", sample)) %>% 
  mutate(total=total/4) %>% slice(1:9)
temporalDS<-read.table(here("data","countreadsDS_ch.tsv"), col.names = c("sample","total"), sep= " ")  %>% mutate(sample = gsub("DS_.*$","final_combined", sample)) %>% mutate(total=total/4)
#temporal_total<-read_tsv(here("data","countreads_ch.tsv.tsv"), col_names = c("sample","total"))
mw_total<-read_tsv(here("data","countreads_mw.tsv.tsv"), col_names = c("sample","total")) 
```

## PREPARE TABLE OF READ COUNTS ASSOCIATED WITH RESISTANCE GENES
```{r}
mw_data <- mw_data %>%
  mutate(
    DS = gsub("^[^_]+_", "", Muestra),
    Muestra = gsub("(\\d+)Gb", "MA_DS\\1_Crowfoot", Muestra)
  )
l_data<-l_data %>% 
  mutate(DS="1405Gb")
t_data<-t_data %>% mutate(DS="5.2Gb",Muestra=gsub("_.*$","_final_combined",Muestra))
last<-last %>% mutate(DS="3.8GB") %>% 
    mutate(Muestra = str_replace(Muestra, "barcode07_DS_2.7GB", "barcode07_DS_3.8GB"))
all<-bind_rows(
  L = l_data, 
  MW = mw_data, 
  temporal=t_data,
  lig_p2=last,
  .id = "folder"
) %>% rename(sample=Muestra)
all_filtered <- all %>%
filter(folder != "MW" | (folder == "MW" & sample %in% sprintf("barcode%02d_MA_DS9_Crowfoot", 1:15))) %>% select(-c(3:5,8:10)) %>%  separate(Feature, into = c("Gen", "Drug_Class"), sep = "_", extra = "merge") %>% 
  separate_rows(Drug_Class, sep = ";") %>% 
   mutate(Drug_Class = str_remove(Drug_Class, "^_+"))
keywords <- c("carbapenem")
all_filtered2<-all_filtered %>% 
  filter(grepl("carbapenem", Drug_Class)) %>% 
  group_by(sample) %>% summarize(conteo=sum(Conteo))
```

## MERGE READ COUNTS WITH LIBRARY SIZES AND METADATA: Barplot
Note: here I used the mean (not the median), and it does not include the seq1212 or seq2111 dataframes
```{r}
table<-bind_rows(mw_total, temporalDS, last_seqDS) %>% inner_join(all_filtered2,by = "sample") %>% mutate(  # Total sum of counts per sample
    fraction = (conteo / total) * 100  # Percentage calculation
  ) %>%
  select(-c(total)
) %>% inner_join(metadata, by="sample") %>% unite("identifier", place, date, sample_id, sep = "_", remove = FALSE) %>% 
    mutate(
    fraction = as.numeric(fraction),
    date = dmy(date)
  ) %>% filter(sample_id != "Mock Community") %>%
  mutate(identifier = gsub("AR{0,1}\\d+(_\\d+)?", "\\1", identifier)) %>%
  mutate(identifier = gsub("_+", "_", identifier)) %>%
  mutate(identifier = gsub("^_|_$", "", identifier)) %>%
  mutate(identifier = gsub("^[[:alnum:]]+_", "", identifier)) %>%
  mutate(
    # Redefine 'identifier' according to your new mapping
    identifier = case_when(
      sample_id == "Mock Community" ~ "Mock Community",
      place == "HMA" & date == as.Date("2024-02-06") ~ "HMA_baseline",
      place == "HMA" & date == as.Date("2024-02-13") ~ "HMA_1_week",
      place == "HMA" & date == as.Date("2024-02-20") ~ "HMA_2_weeks",
      place == "HMA" & date == as.Date("2024-03-12") ~ "HMA_5_weeks",
      place == "HMA" & date == as.Date("2024-05-22") ~ "HMA_3_months",
      place == "INSN" & date %in% as.Date("2024-05-06") ~ "INSN-I_I",
      place == "INSN" & date == as.Date("2024-05-27")~
        "INSN-I_II",
      place == "INSN" & date == as.Date("2024-07-08") ~ "INSN-II",
      TRUE ~ place
    ),
    # Redefine 'group' according to the new 'identifier'
    group = case_when(
      identifier %in% c("INSN-I_I","INSN-I_II", "INSN-II") ~ "INSN",
      identifier == "HCH" ~ "HCH",
      identifier %in% c("HMA_baseline", "HMA_1_week", "HMA_2_weeks", "HMA_5_weeks", "HMA_3_months") ~ "HMA",
      TRUE ~ "Others"
    )
  )

# Order 'identifier' by 'group' and 'date'
table <- table %>%
  arrange(factor(group, levels = c("INSN", "HCH","HMA")), date)

# Update the levels of the 'identifier' factor in the desired order
table$identifier <- factor(table$identifier, levels = unique(table$identifier))

# Calculate 'mean_fraction' and 'mean_directly_reads' by 'identifier'
mean_fraction_df <- table %>%
  group_by(identifier, group) %>%
  summarise(
    mean_fraction = mean(fraction),
    mean_directly_reads = mean(conteo),
    .groups = 'drop'
  )

# Ensure that the levels of the 'identifier' factor match
mean_fraction_df$identifier <- factor(mean_fraction_df$identifier, levels = levels(table$identifier))
esp<-"carbapenem"


# Split the data into two dataframes: one for HMA and one for INSN and HCH
combined_df2_hma <- table %>% filter(group == "HMA")
mean_fraction_df2_hma <- mean_fraction_df %>% filter(group == "HMA")

combined_df2_others <- table %>% filter(group %in% c("INSN", "HCH"))
mean_fraction_df2_others <- mean_fraction_df %>% filter(group %in% c("INSN", "HCH"))
# Define the color palette
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
  "Negative control" = "#CB181D" # Dark red
)

# Ensure that the colors correspond to the levels of 'identifier'

# Create the plot for HMA and store it in 'p_hma'
p_hma <- ggplot() +
  geom_bar(
    data = mean_fraction_df2_hma,
    aes(x = identifier, y = mean_fraction, fill = identifier),
    stat = "identity",
    width = 0.6,
    show.legend = FALSE
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
      c(0, 0.05, 0.1, 0.15, 0.2, seq(1, round(max(combined_df2_hma$fraction))))
    },
    labels = function(breaks) {
      sapply(breaks, function(x) ifelse(x < 0.5, sprintf("%.2f", x), sprintf("%.1f", x)))
    },
    expand = expansion(mult = c(0, 0.10))
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = bquote("Mean percentage of " * italic(.(esp)) * " reads - HMA"),
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
    show.legend = FALSE
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
  # Adjust the Y axis for HCH and INSN
  scale_y_continuous(
    breaks = c(0.05, 0.1, 0.15, 0.2),
    expand = expansion(mult = c(0, 0.10))
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = bquote("Mean percentage of " * italic(.(esp)) * " reads - INSN and HCH"),
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
ggsave("barplot_carbapenem_resistance_hma.png", plot = p_hma, width = 8, height = 6, dpi = 300, bg = "white")
ggsave("barplot_carbapenem_resistance_others.png", plot = p_others, width = 8, height = 6, dpi = 300, bg = "white")


```

## MERGE COUNTS WITH LIBRARY SIZES AND METADATA: Temporal variation line plot
Note: here I did use the median (still missing to add seq1212 or seq2111)
```{r}
keywords <- c("penam", "tetracycline", "fluoroquinolone", "cephalosporin", "aminoglycoside", "carbapenem")
all_filtered2<-all_filtered %>% 
  filter(grepl(paste(keywords, collapse = "|"), Drug_Class)) %>% 
  group_by(sample, Drug_Class) %>% summarize(conteo=sum(Conteo))

table<-bind_rows(mw_total, temporalDS, last_seqDS) %>% inner_join(all_filtered2,by = "sample") %>% mutate(  # Total sum of counts per sample
    fraction = (conteo / total) * 100  # Percentage calculation
  ) %>%
  select(-c(total)
) %>% inner_join(metadata, by="sample")  %>%  filter(
     !str_detect(sample_id, "Mock Community") & 
    !str_detect(sample_id, "Negative Control")
   ) %>%  
    mutate(
    fraction = as.numeric(fraction),
    date = dmy(date)
  )
# Summarize the data for the lines (median fraction by date, place, and species)
data_summary <- table %>%
  group_by(date, place, Drug_Class) %>%
  summarise(
    median_fraction = median(fraction, na.rm = TRUE),  # Median
    q1 = quantile(fraction, 0.25, na.rm = TRUE),  # Lower quartile
    q3 = quantile(fraction, 0.75, na.rm = TRUE),   # Upper quartile
    mean_directly_reads = mean(conteo)
  ) %>%
  ungroup()
gg <- ggplot(table, aes(x = date, y = fraction, color = place, shape = place)) +
  # Median line
  geom_line(
    data = data_summary,
    aes(y = median_fraction, group = place),
    size = 0.8
  ) +
  # Median points with quartiles (Q1 and Q3)
  geom_pointrange(
    data = data_summary,
    aes(y = median_fraction, ymin = q1, ymax = q3),
    position = position_dodge(width = 0.7)
  ) +
   facet_wrap(~Drug_Class, scales = "free_y", labeller = labeller(Drug_Class = function(x) {
    # Remove underscores and convert to title
    str_to_title(gsub("_", " ", x))
  }))+
  labs(
        title = "Percentage abundance of reads assigned to antibiotic resistance genes classes",
    x = "Date",
    y = "Relative abundance (%)",
    color = "Place",
    shape = "Place"
  ) +
  theme_minimal() +
  theme(
     plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank()
  )

ggsave("temporal_variation_resistance.png", plot = gg, width = 8, height = 6, dpi = 300, bg = "white")

```

## ADD CLINICAL SAMPLES TO THE TEMPORAL VARIATION LINE PLOT OF CARBAPENEMS
```{r}
clinical<-read.xlsx(here("data","clinical_HMA_procesados.xlsx"),cols = 1:10, rows = 1:59, colNames = TRUE)
#Fix dates
fecha_base <- as.Date("1900-01-01") - 2  # Subtract 2 days to correct the "1900 error" in Excel
clinical$FECHA_RECEP_MX <- as.Date(clinical$FECHA_RECEP_MX, origin = fecha_base)
clinical$DATE_ANTIBIO <- as.Date(clinical$DATE_ANTIBIO, origin = fecha_base)

#Remove rows with NAs
clinical<-clinical[rowSums(is.na(clinical))==0,]
#Find indices of informative columns
s_i_r_cols <- grep("S/I/R", names(clinical))
#Filter to keep only rows with at least one R in the columns of interest
clinical<-clinical[rowSums(sapply(clinical[s_i_r_cols], str_detect,"R"))>0,]
dim(clinical)
species<-c("Klebsiella pne","Acinetobacter","Pseudomona","Escherichia","Providencia ret","Enterobacter")
clinical_filtered<-clinical %>% select(3,5) %>% filter(str_detect(BACTERIA,str_c(species, collapse = "|")))

table<-bind_rows(mw_total, temporalDS, last_seqDS) %>% inner_join(all_filtered2,by = "sample") %>% mutate(  # Total sum of counts per sample
    fraction = (conteo / total) * 100  # Percentage calculation
  ) %>%
  select(-c(total)
) %>% inner_join(metadata, by="sample")  %>%  filter(
     !str_detect(sample_id, "Mock Community") & 
    !str_detect(sample_id, "Negative Control")
   ) %>%  
    mutate(
    fraction = as.numeric(fraction),
    date = dmy(date)
  )
# Summarize the data for the lines (median fraction by date, place, and species)
data_summary <- table %>%
  group_by(date, place) %>%
  summarise(
    median_fraction = median(fraction, na.rm = TRUE),  # Median
    q1 = quantile(fraction, 0.25, na.rm = TRUE),  # Lower quartile
    q3 = quantile(fraction, 0.75, na.rm = TRUE),   # Upper quartile
    mean_directly_reads = mean(conteo)
  ) %>%
  ungroup()
gg <- ggplot(table, aes(x = date, y = fraction, color = place, shape = place)) +
  # Median line
  geom_line(
    data = data_summary,
    aes(y = median_fraction, group = place),
    size = 0.8
  ) +
  # Median points with quartiles (Q1 and Q3)
  geom_pointrange(
    data = data_summary,
    aes(y = median_fraction, ymin = q1, ymax = q3),
    position = position_dodge(width = 0.7)
  ) +
  labs(
        title = "Percentage abundance of reads assigned to carbapenem resistance genes",
    x = "Date",
    y = "Relative abundance (%)",
    color = "Place",
    shape = "Place"
  ) +
  theme_minimal() +
  theme(
     plot.title = element_text(hjust = 0.5, face = "bold"),
    panel.grid.major = element_line(color = "gray85"),
    panel.grid.minor = element_blank()
  )
min_date <- min(clinical_filtered$FECHA_RECEP_MX)
max_date <- max(clinical_filtered$FECHA_RECEP_MX)
ggc <- gg +
  # Density with limited range and without altering the original scale
  geom_histogram(
    data = clinical_filtered,
    aes(x = FECHA_RECEP_MX,y = after_stat(count * 0.02)),
                 bins = 50,  # adjust the number of bins as needed
                 alpha = 0.5, inherit.aes = FALSE) +
    scale_y_continuous(
    sec.axis = sec_axis(~./0.02, name = "Clinical Resistance Count")
  ) +
  labs(subtitle = "With histogram of resistant clinical bacteria")


ggsave("temporal_variation_resistance+clinical.png", plot = ggc, width = 8, height = 6, dpi = 300, bg = "white")


```




