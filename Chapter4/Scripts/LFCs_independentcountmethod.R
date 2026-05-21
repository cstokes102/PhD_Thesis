#sgcounts LFCs
#FINAL MANUAL LFCS; LL and UL v T0 
rm (list = ls())
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/sgcounts/")
library(dplyr)
library(VennDiagram)
library(tidyverse)
library(ggplot2)
library(stringr)
library(readr)
library(tidyr)
rawsgcounts = read.csv("allsgcounts_withsgREF.csv")
#rename columns so more simple
colnames(rawsgcounts) = c("Guide", "Gene", "T0", "D7", "UL", "LL")
#Convert count columns to numeric
numeric_cols1 <- c("T0", "D7", "UL", "LL")
rawsgcounts[numeric_cols1] <- lapply(rawsgcounts[numeric_cols1], function(x) as.numeric(as.character(x)))
#how many guides exist for each gene in the library
guides_per_gene <- rawsgcounts %>%
  group_by(Gene) %>%
  summarise(total_guides = n(), .groups = "drop")
#go back to counts data and treat anything <10 as 0
rawcounts_strict <- rawsgcounts %>%
  mutate(across(c(T0, D7, UL, LL), ~ ifelse(. < 10, 0, .)))
#look for total guide dropout
totaldropout <- rawcounts_strict %>%
  filter(T0 == 0 & D7 == 0 & UL == 0 & LL == 0)
#419 dropouts
counts_filtered <- anti_join(rawcounts_strict, totaldropout, by = "Guide")
#12578 guides remaining
#add +1 pseudocount
counts_filtered <- counts_filtered %>%
  mutate(across(all_of(numeric_cols1), ~ ifelse(. == 0, 1, .)))
#now compute LFC for each guide; UL v T0
lfc_df <- counts_filtered %>%
  mutate(LFC = log2(UL / T0))
#to account for systematic biases, centre the LFC values using negative control sgRNAs
#Compute the median LFC of negative controls
median_negctrl <- median(
  lfc_df %>% filter(Gene == "negative_control") %>% pull(LFC),
  na.rm = TRUE)
#-1.58
summary(lfc_df$LFC[lfc_df$Gene == "negative_control"])
#currently there are 246 neg control guides
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-12.0351  -9.6357  -1.5800  -4.3740   0.4474   3.0253 
#now summarise LFCs in general
summary(lfc_df$LFC)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-12.4504  -9.7665  -3.8261  -4.5036   0.4562  12.0678
# Subtract the median of controls from all guides
lfc_df <- lfc_df %>%
  mutate(LFC_corrected = LFC - median_negctrl)
summary(lfc_df$LFC_corrected[lfc_df$Gene == "negative_control"])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-10.455  -8.056   0.000  -2.794   2.027   4.605  
#median-centred
#now look at general LFCs now after centring
summary(lfc_df$LFC_corrected)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-10.870  -8.187  -2.246  -2.924   2.036  13.648  
gene_summary <- lfc_df %>%
  group_by(Gene) %>%
  summarise(
    mean_LFC = mean(LFC_corrected, na.rm = TRUE),
    sd_LFC = sd(LFC_corrected, na.rm = TRUE),
    n_guides = n()
  ) %>%
  arrange(desc(mean_LFC))  # optional: sort by effect size
#2314 genes
#look at controls
gene_summary %>%
  filter(Gene == "negative_control")
#Gene             mean_LFC sd_LFC n_guides
#negative_control    -2.79   5.09      246
#quite a big SD for negative control LFCs...
#keep only genes that have at least 3 guides
gene_summary_filtered <- gene_summary %>%
  filter(n_guides >= 3)
#there are now 2294 genes
#use the within-gene variation of guide LFCs to standardize mean effects
z_df <- lfc_df %>%
  filter(Gene %in% gene_summary_filtered$Gene) %>%
  group_by(Gene) %>%
  summarise(
    mean_LFC = mean(LFC, na.rm = TRUE),
    sd_LFC = sd(LFC, na.rm = TRUE),
    n_guides = n()
  ) %>%
  mutate(Zscore = mean_LFC / (sd_LFC / sqrt(n_guides)))
#this is basically a t stat based on within-gene guide variability
#now look for ones with little variation (sd < 2)
UL = z_df %>%
  filter(z_df$sd_LFC < 2)
#now have 179 genes
write.csv(UL, "FINAL_ULhits.csv", row.names = FALSE)
#look at how many of these are just neuron dropouts
neurondropout <- counts_filtered %>% filter(D7 == 1 & UL == 1 & LL == 1)
#2348 guides
neurondropout_summary <- neurondropout %>%
  group_by(Gene) %>%
  summarise(dropout_guides = n(), .groups = "drop")
#1429 genes affected
dropout_vs_total <- guides_per_gene %>%
  left_join(neurondropout_summary, by = "Gene") %>%
  mutate(dropout_guides = ifelse(is.na(dropout_guides), 0, dropout_guides)) %>%
  arrange(desc(dropout_guides))
#how many guides remaining
dropout_vs_total <- guides_per_gene %>%
  left_join(neurondropout_summary, by = "Gene") %>%
  mutate(
    dropout_guides = ifelse(is.na(dropout_guides), 0, dropout_guides),
    guides_remaining = total_guides - dropout_guides
  ) %>%
  arrange(guides_remaining, desc(dropout_guides))
#what genes are no longer represented because of this?
genes_toremove <- dropout_vs_total %>%
  filter(guides_remaining <= 2)
#138 genes with 2 or fewer guides remaining in neurons (not informative/reliable)
#how many of these neuronal dropouts are in my UL hits?
xx = intersect(UL$Gene, genes_toremove$Gene)
#45 overlapping genes
#remove these genes from the hit list
UL_filtered <- UL %>% filter(!Gene %in% genes_toremove$Gene)
#134 genes left
#sort by direction - are they plus or minus LFC compared to T0?
UPos = UL_filtered %>% filter(mean_LFC > 0) %>% pull(Gene)
#61 positives - enriched relative to T0
UNeg = UL_filtered %>% filter(mean_LFC < 0) %>% pull(Gene)
#73 negatives - depleted relative to T0
write.csv(data.frame(UPos), "FINAL_ULhitspos.csv", row.names = FALSE)
write.csv(data.frame(UNeg), "FINAL_ULhitsneg.csv", row.names = FALSE)

#now do same for PSL v T0 LFCs
#go through the same up until this:
#now compute LFC for each guide; LL v T0
lfc_df <- counts_filtered %>%
  mutate(LFC = log2(LL / T0))
#to account for systematic biases, centre the LFC values using negative control sgRNAs
#Compute the median LFC of negative controls
median_negctrl <- median(
  lfc_df %>% filter(Gene == "negative_control") %>% pull(LFC),
  na.rm = TRUE)
# -9.288862
summary(lfc_df$LFC[lfc_df$Gene == "negative_control"])
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-11.7398 -10.2891  -9.2889  -6.2755   0.4904   4.1983 
# Subtract the median of controls from all guides
lfc_df <- lfc_df %>%
  mutate(LFC_corrected = LFC - median_negctrl)
summary(lfc_df$LFC_corrected[lfc_df$Gene == "negative_control"])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-2.451  -1.000   0.000   3.013   9.779  13.487 
#median-centred
gene_summary <- lfc_df %>%
  group_by(Gene) %>%
  summarise(
    mean_LFC = mean(LFC_corrected, na.rm = TRUE),
    sd_LFC = sd(LFC_corrected, na.rm = TRUE),
    n_guides = n()
  ) %>%
  arrange(desc(mean_LFC))  # optional: sort by effect size
#2314 genes
#look at controls
gene_summary %>%
  filter(Gene == "negative_control")
#Gene             mean_LFC sd_LFC n_guides
#negative_control     3.01   5.34      246
#quite a big SD for negative control LFCs...
#keep only genes that have at least 3 guides
gene_summary_filtered <- gene_summary %>%
  filter(n_guides >= 3)
#unsure where it is getting the number of guides from - assuming from T0?
#there are now 2294 genes
#use the within-gene variation of guide LFCs to standardize mean effects
z_df <- lfc_df %>%
  filter(Gene %in% gene_summary_filtered$Gene) %>%
  group_by(Gene) %>%
  summarise(
    mean_LFC = mean(LFC, na.rm = TRUE),
    sd_LFC = sd(LFC, na.rm = TRUE),
    n_guides = n()
  ) %>%
  mutate(Zscore = mean_LFC / (sd_LFC / sqrt(n_guides)))
#this is basically a t stat based on within-gene guide variability
#now look for ones with little variation (sd < 2)
LL = z_df %>%
  filter(z_df$sd_LFC < 2)
#now have 367 genes
#look at how many of these are just neuron dropouts
neurondropout <- counts_filtered %>% filter(D7 == 1 & UL == 1 & LL == 1)
#2348 guides
neurondropout_summary <- neurondropout %>%
  group_by(Gene) %>%
  summarise(dropout_guides = n(), .groups = "drop")
#1429 genes affected
dropout_vs_total <- guides_per_gene %>%
  left_join(neurondropout_summary, by = "Gene") %>%
  mutate(dropout_guides = ifelse(is.na(dropout_guides), 0, dropout_guides)) %>%
  arrange(desc(dropout_guides))
#how many guides remaining
dropout_vs_total <- guides_per_gene %>%
  left_join(neurondropout_summary, by = "Gene") %>%
  mutate(
    dropout_guides = ifelse(is.na(dropout_guides), 0, dropout_guides),
    guides_remaining = total_guides - dropout_guides
  ) %>%
  arrange(guides_remaining, desc(dropout_guides))
#what genes are no longer represented because of this?
genes_toremove <- dropout_vs_total %>%
  filter(guides_remaining <= 2)
#138 genes with 2 or fewer guides remaining in neurons (not informative/reliable)
#how many of these neuronal dropouts are in my LL hits?
aa = intersect(LL$Gene, genes_toremove$Gene)
#68 overlapping genes
#remove these genes from the hit list
LL_filtered <- LL %>% filter(!Gene %in% genes_toremove$Gene)
#299 genes left
write.csv(LL_filtered, "FINAL_LLhits.csv", row.names = FALSE)
#sort by direction - are they plus or minus LFC compared to T0?
LPos = LL_filtered %>% filter(mean_LFC > 0) %>% pull(Gene)
#14 positives - enriched relative to T0
LNeg = LL_filtered %>% filter(mean_LFC < 0) %>% pull(Gene)
#285 negatives - depleted relative to T0
write.csv(data.frame(LPos), "FINAL_LLhitspos.csv", row.names = FALSE)
write.csv(data.frame(LNeg), "FINAL_LLhitsneg.csv", row.names = FALSE)