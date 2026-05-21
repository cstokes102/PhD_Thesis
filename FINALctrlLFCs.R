rm (list = ls())
library(dplyr)
library(VennDiagram)
library(tidyverse)
library(ggplot2)
library(stringr)
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/all_counts_controlnorm/")
rawcounts = read.delim("all_controlnorm.count_normalized.txt")
#rename columns so more simple
colnames(rawcounts) = c("Guide", "Gene", "T0", "D7", "UL", "LL")
#Convert count columns to numeric
numeric_cols1 <- c("T0", "D7", "UL", "LL")
rawcounts[numeric_cols1] <- lapply(rawcounts[numeric_cols1], function(x) as.numeric(as.character(x)))
#how many guides exist for each gene in the library
guides_per_gene <- rawcounts %>%
  group_by(Gene) %>%
  summarise(total_guides = n(), .groups = "drop")
#go back to counts data and treat anything <10 as 0
rawcounts_strict <- rawcounts %>%
  mutate(across(c(T0, D7, UL, LL), ~ ifelse(. < 10, 0, .)))
#look for total guide dropout
totaldropout <- rawcounts_strict %>%
  filter(T0 == 0 & D7 == 0 & UL == 0 & LL == 0)
#502 dropouts
counts_filtered <- anti_join(rawcounts_strict, totaldropout, by = "Guide")
#12506 guides remaining
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
#2.525533
# Subtract the median of controls from all guides
lfc_df <- lfc_df %>%
  mutate(LFC_corrected = LFC - median_negctrl)
summary(lfc_df$LFC_corrected[lfc_df$Gene == "negative_control"])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-11.567  -9.092   0.000  -3.322   1.900   4.529   
#median-centred
#now look at general LFCs now after centring
summary(lfc_df$LFC_corrected)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-11.983  -9.259  -2.684  -3.532   1.925  10.581   
gene_summary <- lfc_df %>%
  group_by(Gene) %>%
  summarise(
    mean_LFC = mean(LFC, na.rm = TRUE),
    sd_LFC = sd(LFC, na.rm = TRUE),
    n_guides = n()
  ) %>%
  arrange(desc(mean_LFC))  # optional: sort by effect size
#2313 genes
#look at controls
gene_summary %>%
  filter(Gene == "negative_control")
#Gene             mean_LFC sd_LFC n_guides
#negative_control   -3.32   5.55      246
#quite a big SD for negative control LFCs...
#keep only genes that have at least 3 guides
gene_summary_filtered <- gene_summary %>%
  filter(n_guides >= 3)
#there are now 2287 genes
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
#now have 174 genes
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
#86 genes with 2 or fewer guides remaining in neurons (not informative/reliable)
#how many of these neuronal dropouts are in my UL hits?
xx = intersect(UL$Gene, genes_toremove$Gene)
#32 overlapping genes
#remove these genes from the hit list
UL_filtered <- UL %>% filter(!Gene %in% genes_toremove$Gene)
#142 genes left
write.csv(UL_filtered, "FINAL_ULhits.csv", row.names = FALSE)
#sort by direction - are they plus or minus LFC compared to T0?
UPos = UL_filtered %>% filter(mean_LFC > 0) %>% pull(Gene)
#71 positives - enriched relative to T0
UNeg = UL_filtered %>% filter(mean_LFC < 0) %>% pull(Gene)
#71 negatives - depleted relative to T0
write.csv(data.frame(UPos), "FINAL_ULhitspos.csv", row.names = FALSE)
write.csv(data.frame(UNeg), "FINAL_ULhitsneg.csv", row.names = FALSE)

#now do same for LL v T0 LFCs
#look for total guide dropout, but dont do the <10 reads threshold when this is normalised already
totaldropout <- rawcounts %>%
  filter(T0 == 0 & D7 == 0 & UL == 0 & LL == 0)
#433
counts_filtered <- anti_join(rawcounts, totaldropout, by = "Guide")
#12506 guides remaining
#add +1 pseudocount
counts_filtered <- counts_filtered %>%
  mutate(across(all_of(numeric_cols1), ~ ifelse(. == 0, 1, .)))
#now compute LFC for each guide; LL v T0
lfc_df <- counts_filtered %>%
  mutate(LFC = log2(LL / T0))
# #to account for systematic biases, centre the LFC values using negative control sgRNAs
# #Compute the median LFC of negative controls
# median_negctrl <- median(
#   lfc_df %>% filter(Gene == "negative_control") %>% pull(LFC),
#   na.rm = TRUE)
# # -0.4435
# # Subtract the median of controls from all guides
# lfc_df <- lfc_df %>%
#   mutate(LFC_corrected = LFC - median_negctrl)
#median-centred
gene_summary <- lfc_df %>%
  group_by(Gene) %>%
  summarise(
    mean_LFC = mean(LFC, na.rm = TRUE),
    sd_LFC = sd(LFC, na.rm = TRUE),
    n_guides = n()
  ) %>%
  arrange(desc(mean_LFC))  # optional: sort by effect size
#2313 genes
#look at controls
# gene_summary %>%
#   filter(Gene == "negative_control")
#Gene             mean_LFC sd_LFC n_guides
#negative_control     1.59   6.48      246
# big SD for negative control LFCs...
#keep only genes that have at least 3 guides
gene_summary_filtered <- gene_summary %>%
  filter(n_guides >= 3)
#unsure where it is getting the number of guides from - assuming from T0?
#there are now 2287 genes
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
#now have 85 genes
#look at how many of these are just neuron dropouts
neurondropout <- counts_filtered %>% filter(D7 == 1 & UL == 1 & LL == 1)
neurondropout_summary <- neurondropout %>%
  group_by(Gene) %>%
  summarise(dropout_guides = n(), .groups = "drop")
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
#86 genes with 2 or fewer guides remaining in neurons (not informative/reliable)
#how many of these neuronal dropouts are in my LL hits?
aa = intersect(LL$Gene, genes_toremove$Gene)
#17 overlapping genes
#remove these genes from the hit list
LL_filtered <- LL %>% filter(!Gene %in% genes_toremove$Gene)
#67 genes left
write.csv(LL_filtered, "FINAL_LLhits.csv", row.names = FALSE)
#sort by direction - are they plus or minus LFC compared to T0?
LPos = LL_filtered %>% filter(mean_LFC > 0) %>% pull(Gene)
#18 positives - enriched relative to T0
LNeg = LL_filtered %>% filter(mean_LFC < 0) %>% pull(Gene)
#49 negatives - depleted relative to T0
write.csv(data.frame(LPos), "FINAL_LLhitspos.csv", row.names = FALSE)
write.csv(data.frame(LNeg), "FINAL_LLhitsneg.csv", row.names = FALSE)
#why is it behaving weird?
ggplot(lfc_df %>% filter(Gene == "negative_control"), aes(LFC)) +
  geom_histogram(bins = 50) +
  theme_classic()
#there seems to be 3 clusters of negative control LFCs
#should not manually median-centre these again