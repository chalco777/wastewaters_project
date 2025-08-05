---
title: "Resistome Carbapenem Analysis"
output: html_document
date: "2024-11-05"
---

## What is this notebook about?

This notebook aggregates carbapenem-gene counts across runs, converts to CPM using sample-specific totals, applies log2(CPM+1), computes Bray–Curtis distances, and performs PCoA to ordain samples, coloring points by biologically meaningful identifiers derived from place+date (e.g., HMA_baseline, INSN-I/II, HCH, Mock). 
It then ranks drivers of separation by correlating each gene’s CPM profile with a PCoA axis (Spearman, BH-FDR) and exporting the top 10 genes. 
Other analysis include a batch-colored PCoA and an exploratory category-level PCoA using Drug_Class abundances.
I built the CPM matrix, ran the ordination, curated palettes/labels for interpretability, generated the top-genes table, and added batch and category analyses.
I also included a raw ConQuR pass to illustrate batch-effect correction on the carbapenem gene table.

## Import libraries as usually

```{r}
library(tidyverse)
library(vegan)
library(glue)
library(ggrepel)
library(openxlsx)
library(doParallel)
library(ConQuR)
library(ggtext)
setwd("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/loui/pcoa_resistome")
mw_data<-read_tsv("count_mw_allDS.tsv")
l_data<-read_tsv("count_L_1405DS.tsv")
t_data<-read_tsv("count_wastewatersch_5.2DS.tsv")
last<-read_tsv("count_LIG_P2_3.8DS.tsv")
```

## PREPARE METADATA

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
metadata<-bind_rows(metadata,last_temp) %>% 
  mutate(place = case_when(sample_id=="Mock Community" ~ "Mock Community",str_detect(sample_id, "Negative")~ "Negative Control", TRUE ~ place)) %>% 
  filter(!str_detect(sample_id, "Negative"))
```

## PREPARE COUNTS

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
all_filtered2<-all_filtered %>% 
  filter(Drug_Class=="carbapenem") %>% 
  group_by(sample,Gen) %>% summarize(conteo=sum(Conteo))
```

## FIND CPM AND CREATE THE PCOA

```{r}
###ADDING TOTAL COUNTS TO FIND CPM
#last_seq_total<-read_tsv("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/loui/abundance/countreads_17sset.tsv.tsv", col_names = c("sample","total"))%>% mutate(sample = paste0(sample,"_DS_3.8GB"))

last_seqDS<-read.table("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/loui/abundance/countreadsDS_17sset.tsv.tsv", col.names = c("sample","total"), sep= " ")%>% mutate(sample = gsub(".fastq.gz.*$","", sample)) %>% 
  mutate(total=total/4) %>% slice(1:9)

temporalDS<-read.table("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/loui/abundance/countreadsDS_ch.tsv.txt", col.names = c("sample","total"), sep= " ")  %>% mutate(sample = gsub("DS_.*$","final_combined", sample)) %>% mutate(total=total/4)

#temporal_total<-read_tsv("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/loui/abundance/countreads_ch.tsv.tsv", col_names = c("sample","total"))
mw_total<-read_tsv("C:/Users/DAVID 21/OneDrive/Documentos/Mirkoslab/loui/abundance/countreads_mw.tsv.tsv", col_names = c("sample","total")) 


cpm_norm<-bind_rows(mw_total, temporalDS, last_seqDS) %>% inner_join(all_filtered2,by = "sample") %>% mutate(  # Total sum of counts per sample
    cpm = (conteo / total) * 1e6  # CPM calculation
  ) %>%
  select(-c(total,conteo)
)
cpm_matrix<-cpm_norm %>% pivot_wider(names_from = Gen, values_from = cpm, values_fill = 0) %>%
  column_to_rownames("sample") %>%
  as.matrix()
log_cpm_matrix<-log2(cpm_matrix+1)
#sqr_cpm_matrix<-sqrt(cpm_matrix)
# Calculate Bray-Curtis distance
#bray_dist <- vegdist(cpm_matrix, method = "bray")
bray_dist <- vegdist(log_cpm_matrix, method = "bray")

# Perform PCoA
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2, add=TRUE)
positions <- pcoa$points #extract the vectors
colnames(positions) <- c("pcoa1", "pcoa2") #assign names to coordinates
positions<-positions[order(rownames(positions)), ] #order by row name
#find explained percentage
expl<-(pcoa$eig / sum(pcoa$eig))*100
exp<-format(round(expl[1:2],digits=1),nsmall=1, trim=TRUE)
labs<-c(glue("PCo 1 ({exp[1]}%)"),
        glue("PCo 2 ({exp[2]}%)"))
# Convert the ordered PCoA results to tibble and join with df_variable
# Improved PCoA plot with sample labels
palette_colors <- c(
  "HMA_baseline"   = "#2171B5",  # Darker blue
  "HMA_1_week"     = "#4292C6",  # Medium blue
  "HMA_2_weeks"    = "#6BAED6",  # Lighter blue
  "HMA_5_weeks"    = "#9ECAE1",  # Even lighter blue
  "HMA_3_months"   = "#2131B5",  
  "INSN-I"           = "#74C476",
  "INSN-II" = "#A1D99B",
    "HCH"            = "#6A51A3",  # Dark purple
  "Mock Community" = "#F16913",  # Medium orange
  "Negative control" = "#CB181D" # Dark red
)

pc<-positions %>%
  as_tibble(rownames = "sample") %>%
  inner_join(metadata, by = "sample") %>% 
  #filter(!is.na(place)) %>%
  #filter(sample_id!="Mock Community") %>% 
  mutate(
    date = dmy(date)
  ) %>%
  mutate(replicate=gsub("^[^_]+_","",sample_id)) %>% 
  mutate(
    identifier = case_when(sample_id=="Mock Community" ~ "Mock Community",
      place == "HMA" & date == as.Date("2024-02-06") ~ "HMA_baseline",
      place == "HMA" & date == as.Date("2024-02-13") ~ "HMA_1_week",
      place == "HMA" & date == as.Date("2024-02-20") ~ "HMA_2_weeks",
      place == "HMA" & date == as.Date("2024-03-12") ~ "HMA_5_weeks",
      place == "HMA" & date == as.Date("2024-05-22") ~ "HMA_3_months",
      place == "INSN" & date %in% as.Date(c("2024-05-06", "2024-05-27")) ~ "INSN-I",
      place == "INSN" & date == as.Date("2024-07-08") ~ "INSN-II",
      TRUE ~ place
    )) 

# Create a vector of 'identifier' levels ordered by date
ordered_levels <- pc %>%
  select(identifier, date, place) %>%
  distinct() %>%
  mutate(place = factor(place, levels = c("HMA", "INSN", "HCH","Mock Community", "Negative Control"))) %>% # Specific order of 'place'
  arrange(place, date) %>%
  pull(identifier)
ordered_levels <- unique(ordered_levels)
pc <- pc %>%
  mutate(identifier = factor(identifier, levels = ordered_levels))

gg<-pc%>%ggplot(aes(x = pcoa1, y = pcoa2, color = identifier)) +
  geom_point(size = 3, show.legend = TRUE)+
  #stat_ellipse(aes(fill = identifier), geom = "polygon", alpha = 0.2, level = 0.95, show.legend = FALSE)+
  #stat_ellipse(aes(fill = identifier), geom="polygon",alpha=0.2, show.legend = FALSE)+
  #geom_jitter(size = 2.5, stroke = 0.8,width = 0.005, height = 0.005) +  # Place geom_point without size in aes()
  # geom_text_repel(aes(label = identifier), size = 3, max.overlaps = 2, force=2)+#tras
  labs(color = "Sample", x = labs[1], y = labs[2]) +  
  theme_minimal(base_size = 16) +  
  theme(
    axis.title = element_text(size = 16, face = "bold"),  
    axis.text = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  
  ) +
  scale_color_manual(values = palette_colors) +  
guides(fill = "none")+
  ggtitle("PCoA of Carbapenem Genes")+
  theme(plot.title = element_markdown())


#+ xlim(c(-0.105, -0.095))

```

## Which genes contribute to the variation on our PCoA axes?

```{r}
library(broom)

# Extract PCo1 scores and convert to a dataframe
positions_df <- positions %>%
  as_tibble(rownames = "sample") %>%
  select(sample, pcoa2)
# Convert CPM matrix to a dataframe and reset row names
cpm_df <- as.data.frame(cpm_matrix) %>%
  rownames_to_column("sample")

# Merge PCo1 scores with CPM data
merged_df <- left_join(positions_df, cpm_df, by = "sample")

# Get list of gene names
genes <- colnames(merged_df)[-(1:2)]  # Exclude 'sample' and 'pcoa1' columns

# Compute Spearman correlation between each gene and PCo1 scores
cor_results <- map_dfr(genes, function(gene) {
  cor_test <- cor.test(merged_df[[gene]], merged_df$pcoa2, method = "spearman")
  tibble(
    Gene = gene,
    Correlation = cor_test$estimate,
    P_value = cor_test$p.value
  )
})

# Adjust p-values for multiple testing using Benjamini-Hochberg method
cor_results <- cor_results %>%
  mutate(Adjusted_P_value = p.adjust(P_value, method = "BH"))

# Select the top 10 genes with the highest absolute correlation with PCo1
top_genes <- cor_results %>%
  arrange(desc(abs(Correlation))) %>%
  slice(1:10)
library(ggpubr)

# Format the table for publication
top_genes_formatted <- top_genes %>%
  mutate(
    Correlation = round(Correlation, 3),
    P_value = formatC(P_value, format = "e", digits = 2),
    Adjusted_P_value = formatC(Adjusted_P_value, format = "e", digits = 2)
  )

# Display the top 10 genes in a beautiful table
library(knitr)
library(kableExtra)
# Save the table as HTML
save_kable(
  kable(top_genes_formatted, align = "lccc", caption = "Top 10 Genes Influencing PCo2 in Carbapenem Genes PCoA (log2)") %>%
    kable_styling(full_width = FALSE, font_size = 12),
  file = "top_genes_log2_table_carbapenemgenespcoa.html"
)
#ggsave("top_genes_table.png", tabla_gg, width = 8, height = 4, dpi = 700)

```

## Let's hope there is no grouping by BATCH

```{r}
pc_batch<-pc %>% ggplot(aes(x = pcoa1, y = pcoa2, color = folder)) +
  geom_point(size = 3, show.legend = TRUE)+
  #geom_jitter(size = 2.5, stroke = 0.8,width = 0.005, height = 0.005) +  # Place geom_point without size in aes()
  # geom_text_repel(aes(label = identifier), size = 3, max.overlaps = 2, force=2)+#tras
  labs(color = "Batch", x = labs[1], y = labs[2]) +  
  theme_minimal(base_size = 16) +  
  theme(
    axis.title = element_text(size = 16, face = "bold"),  
    axis.text = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  
  ) +
  scale_color_brewer(palette = "Set1") +  
  ggtitle("PCoA of Carbapenem Genes by Batch")
```

## Now we will build a PCoA from CATEGORIES of Antibiotic Resistance Genes

```{r}
cat_all_filtered2<-all_filtered %>% 
  #filter(Drug_Class=="carbapenem") %>% 
  group_by(sample,Drug_Class) %>% summarize(conteo=sum(Conteo))

cat_cpm_norm<-bind_rows(mw_total, temporal_total, last_seq_total) %>% inner_join(cat_all_filtered2,by = "sample") %>% mutate(  # Total sum of counts per sample
    cpm = (conteo / total) * 1e6  # CPM calculation
  ) %>%
  select(-c(total,conteo)
)
cat_cpm_matrix<-cat_cpm_norm %>% pivot_wider(names_from = Drug_Class, values_from = cpm, values_fill = 0) %>%
  column_to_rownames("sample") %>%
  as.matrix()

# Calculate Bray-Curtis distance
bray_dist <- vegdist(cat_cpm_matrix, method = "bray")

# Perform PCoA
pcoa <- cmdscale(bray_dist, eig = TRUE, k = 2, add=TRUE)


positions <- pcoa$points #extract the vectors
colnames(positions) <- c("pcoa1", "pcoa2") #assign names to coordinates
positions<-positions[order(rownames(positions)), ] #order by row name

#find explained percentage
expl<-(pcoa$eig / sum(pcoa$eig))*100
exp<-format(round(expl[1:2],digits=1),nsmall=1, trim=TRUE)
labs<-c(glue("PCo 1 ({exp[1]}%)"),
        glue("PCo 2 ({exp[2]}%)"))
         
# Convert the ordered PCoA results to tibble and join with df_variable
pc<-positions %>%
  as_tibble(rownames = "sample") %>%
  inner_join(metadata, by = "sample") %>% 
  #filter(!is.na(place)) %>%
  #filter(sample_id!="Mock Community") %>% 
  mutate(
    date = dmy(date)
  ) %>%
  mutate(replicate=gsub("^[^_]+_","",sample_id)) %>% 
  mutate(
    identifier = case_when(
            sample_id=="Mock Community" ~ "Mock Community",
      place == "HMA" & date == as.Date("2024-02-06") ~ "HMA_baseline",
      place == "HMA" & date == as.Date("2024-02-13") ~ "HMA_1_week",
      place == "HMA" & date == as.Date("2024-02-20") ~ "HMA_2_weeks",
      place == "HMA" & date == as.Date("2024-03-12") ~ "HMA_5_weeks",
      place == "HMA" & date == as.Date("2024-05-22") ~ "HMA_3_months",
      TRUE ~ place
    )
  ) %>% ggplot(aes(x = pcoa1, y = pcoa2, color = identifier)) +
  geom_point(size = 2.5, show.legend = TRUE)+
  #stat_ellipse(aes(fill = identifier), geom = "polygon", alpha = 0.2, level = 0.95, show.legend = FALSE)+
  #stat_ellipse(aes(fill = identifier), geom="polygon",alpha=0.2, show.legend = FALSE)+
  #geom_jitter(size = 2.5, stroke = 0.8,width = 0.005, height = 0.005) +  # Place geom_point without size in aes()
  # geom_text_repel(aes(label = identifier), size = 3, max.overlaps = 2, force=2)+#tras
  labs(color = "Sample", x = labs[1], y = labs[2]) +  
  theme_minimal(base_size = 16) +  
  theme(
    axis.title = element_text(size = 16, face = "bold"),  
    axis.text = element_text(size = 14),  
    legend.title = element_text(size = 14),  
    legend.text = element_text(size = 12),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  
  ) +
  scale_color_brewer(palette = "Set1") +  
guides(fill = "none") +  
  ggtitle("PCoA of Categories")
#+ xlim(c(-0.105, -0.095))


```

## Which are the top categories influencing this PCOA?

```{r}

# Extract PCo1 scores and convert to a dataframe
positions_df <- positions %>%
  as_tibble(rownames = "sample") %>%
  select(sample, pcoa1)
# Convert CPM matrix to a dataframe and reset row names
cpm_df <- as.data.frame(cat_cpm_matrix) %>%
  rownames_to_column("sample")

# Merge PCo1 scores with CPM data
merged_df <- left_join(positions_df, cpm_df, by = "sample")

# Calculate the Spearman correlation with the transformed data
cor_results <- map_dfr(genes, function(gene) {
  cor_test <- cor.test(merged_df[[gene]], merged_df_transformed$pcoa1, method = "kendall")
  tibble(
    Gene = gene,
    Correlation = cor_test$estimate,
    P_value = cor_test$p.value
  )
})

# Adjust p-values for multiple testing using Benjamini-Hochberg method
cor_results <- cor_results %>%
  mutate(Adjusted_P_value = p.adjust(P_value, method = "BH"))

# Select the top 10 genes with the highest absolute correlation with PCo1
top_genes <- cor_results %>%
  arrange(desc(abs(Correlation))) %>%
  slice(1:10)
library(ggpubr)

# Format the table for publication
top_genes_formatted <- top_genes %>%
  mutate(
    Correlation = round(Correlation, 3),
    P_value = formatC(P_value, format = "e", digits = 2),
    Adjusted_P_value = formatC(Adjusted_P_value, format = "e", digits = 2)
  )

# Display the top 10 genes in a beautiful table
library(knitr)
library(kableExtra)
# Save the table as HTML
save_kable(
  kable(top_genes_formatted, align = "lccc", caption = "Top 10 Categories Influencing PCo1 in Categories PCoA") %>%
    kable_styling(full_width = FALSE, font_size = 12),
  file = "top_genes_table_categoriespcoaoff.html"
)
#ggsave("top_genes_table.png", tabla_gg, width = 8, height = 4, dpi = 700)
```

## USING RAW CONQUR

```{r}
metadata<-metadata %>% 
    filter(!is.na(place)) %>%
  filter(sample_id!="Mock Community")
my_sample_data= inner_join(all_filtered2, metadata, by="sample")
my_sample_data2<-my_sample_data %>% 
  pivot_wider(names_from = Gen, values_from = conteo, values_fill = 0)
x<-my_sample_data2$sample

my_sample_data2<-my_sample_data2[,-c(1,2)]%>% 
  mutate(across(1:4, as.factor))
rownames(my_sample_data2)<-x

taxa = as.matrix(my_sample_data2[, -(1:4)] )
batchid = my_sample_data2[, 'folder']
summary(batchid)

# Extract the batchids vector
batchid_vector <- batchid$folder

# Create a logical matrix indicating the presence of each taxon (count greater than zero)
present_matrix <- taxa > 0

# Calculate the number of unique batchids in which each taxon is present
num_batchids_per_taxon <- apply(present_matrix, 2, function(present_vec) {
  batchids_present <- batchid_vector[present_vec]
  length(unique(batchids_present))
})

# Select taxa that are present in at least two batchids
taxa_selected <- taxa[, num_batchids_per_taxon >= 2]

# Let's see how many taxa meet the condition
length(num_batchids_per_taxon[num_batchids_per_taxon >= 2])
rownames(taxa_selected)<-x
rownames(covar)<-x
rownames(batchid)<-x
covar = my_sample_data2[, c('place')]
summary(covar)

covar<-covar %>% mutate(place=as.factor(place))


taxa_corrected_raw1=ConQuR(
  tax_tab = taxa_selected,
  batchid = batchid$folder,
  covariates = covar,
  batch_ref = "MW_P2",
  logistic_lasso = FALSE,
  quantile_type = "standard",
  simple_match = FALSE,
  lambda_quantile = "2p/n",
  interplt = FALSE,
  delta = 0.4999,
  taus = seq(0.005, 0.995, by = 0.005),
  num_core = 2
)

```











 
