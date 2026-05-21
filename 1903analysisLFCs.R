#FINAL MANUAL LFCS; LL and UL v T0 
rm (list = ls())
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/")
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(VennDiagram)
#but this time include at the start that any guide with less than 10 reads should be considered as zero
#read in raw counts data
rawcounts = read.delim("allnonorm.count.txt")
#Convert count columns to numeric
numeric_cols <- c("T0", "D7", "UL", "LL")
rawcounts[numeric_cols] <- lapply(rawcounts[numeric_cols], function(x) as.numeric(as.character(x)))
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
#449 dropouts
counts_filtered <- anti_join(rawcounts_strict, totaldropout, by = "sgRNA")
#12559 guides remaining
#add +1 pseudocount
counts_filtered <- counts_filtered %>%
  mutate(across(all_of(numeric_cols), ~ ifelse(. == 0, 1, .)))
#now compute LFC for each guide; UL v T0
lfc_df <- counts_filtered %>%
  mutate(LFC = log2(UL / T0))
#to account for systematic biases, centre the LFC values using negative control sgRNAs
#Compute the median LFC of negative controls
median_negctrl <- median(
  lfc_df %>% filter(Gene == "negative_control") %>% pull(LFC),
  na.rm = TRUE)
#-1.277525
summary(lfc_df$LFC[lfc_df$Gene == "negative_control"])
#currently there are 246 neg control guides
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-11.9226  -9.5274  -1.2775  -4.2141   0.6335   3.2512 
#now summarise LFCs in general
summary(lfc_df$LFC)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-12.3393  -9.6795  -5.3219  -4.3973   0.6443  12.1680 
# Subtract the median of controls from all guides
lfc_df <- lfc_df %>%
  mutate(LFC_corrected = LFC - median_negctrl)
summary(lfc_df$LFC_corrected[lfc_df$Gene == "negative_control"])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-10.645  -8.250   0.000  -2.937   1.911   4.529 
#median-centred
#now look at general LFCs now after centring
summary(lfc_df$LFC_corrected)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-11.062  -8.402  -4.044  -3.120   1.922  13.446 
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
#negative_control    -2.94   5.14      246
#Gene             mean_LFC sd_LFC n_guides
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
UL = z_df %>%
  filter(z_df$sd_LFC < 2)
#now have 176 genes
#sort by direction - are they plus or minus LFC compared to T0?
UPos = UL %>% filter(mean_LFC > 0) %>% pull(Gene)
#66 positives - enriched relative to T0
UNeg = UL %>% filter(mean_LFC < 0) %>% pull(Gene)
#110 negatives - depleted relative to T0
#look at how many of these are just neuron dropouts
neurondropout <- counts_filtered %>% filter(D7 == 0 & UL == 0 & LL == 0)
#2388 guides
neurondropout_summary <- neurondropout %>%
  group_by(Gene) %>%
  summarise(dropout_guides = n(), .groups = "drop")
#1445 genes affected
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
#143 genes with 2 or fewer guides remaining in neurons (not informative/reliable)
#how many of these neuronal dropouts are in my UL hits?
xx = intersect(UL$Gene, genes_toremove$Gene)
#45 overlapping genes
#remove these genes from the hit list
UL_filtered <- UL %>% filter(!Gene %in% genes_toremove$Gene)
#131 genes left
write.csv(UL_filtered, "FINAL_allnonormhits.csv", row.names = FALSE)
#sort by direction - are they plus or minus LFC compared to T0?
UPos = UL_filtered %>% filter(mean_LFC > 0) %>% pull(Gene)
#66 positives - enriched relative to T0
UNeg = UL_filtered %>% filter(mean_LFC < 0) %>% pull(Gene)
#65 negatives - depleted relative to T0
#this means all of those uninformative genes missing in neurons were in the depletion list as expected
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
# -9.19962
summary(lfc_df$LFC[lfc_df$Gene == "negative_control"])
#currently there are 246 neg control guides
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-11.6416 -10.1935  -9.1996  -6.2527   0.4713   4.1861 
#I know 44 of these 246 ctrl guides are missing in the neuronal pops
#unsure if I should get rid so they stop warping things
lfcnegs = lfc_df %>% filter(Gene == "negative_control")
xx = intersect(neurondropout$sgRNA, lfcnegs$sgRNA)

#now summarise LFCs in general
summary(lfc_df$LFC)
#Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#-12.285 -10.243  -9.077  -6.131   0.495  12.623 
# Subtract the median of controls from all guides
lfc_df <- lfc_df %>%
  mutate(LFC_corrected = LFC - median_negctrl)
summary(lfc_df$LFC_corrected[lfc_df$Gene == "negative_control"])
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-2.4420 -0.9939  0.0000  2.9469  9.6709 13.3857 
#median-centred
#now look at general LFCs now after centring
summary(lfc_df$LFC_corrected)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-3.0858 -1.0436  0.1228  3.0690  9.6946 21.8226  
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
#negative_control     2.95   5.28      246
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
    mean_LFC = mean(LFC_corrected, na.rm = TRUE),
    sd_LFC = sd(LFC_corrected, na.rm = TRUE),
    n_guides = n()
  ) %>%
  mutate(Zscore = mean_LFC / (sd_LFC / sqrt(n_guides)))
#this is basically a t stat based on within-gene guide variability
#now look for ones with little variation (sd < 2)
LL = z_df %>%
  filter(z_df$sd_LFC < 2)
#now have 379 genes

#why is it behaving weird?
ggplot(lfc_df %>% filter(Gene == "negative_control"), aes(LFC_corrected)) +
  geom_histogram(bins = 50) +
  theme_classic()
#two peaks
neg_ctrl <- rawcounts %>%
  filter(Gene == "negative_control") %>%
  mutate(LFC = log2((LL + 1)/(T0 + 1)))





#sort by direction - are they plus or minus LFC compared to T0?
LPos = LL %>% filter(mean_LFC > 0) %>% pull(Gene)
#14 positives - enriched relative to T0
LNeg = LL %>% filter(mean_LFC < 0) %>% pull(Gene)
#365 negatives - depleted relative to T0
#look at how many of these are just neuron dropouts
neurondropout <- counts_filtered %>% filter(D7 == 1 & UL == 1 & LL == 1)
#2388 guides
neurondropout_summary <- neurondropout %>%
  group_by(Gene) %>%
  summarise(dropout_guides = n(), .groups = "drop")
#1445 genes affected
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
#143 genes with 2 or fewer guides remaining in neurons (not informative/reliable)
#how many of these neuronal dropouts are in my LL hits?
aa = intersect(LL$Gene, genes_toremove$Gene)
#75 overlapping genes
#remove these genes from the hit list
LL_filtered <- LL %>% filter(!Gene %in% genes_toremove$Gene)
write.csv(LL_filtered, "FINAL_allnonormLLhits.csv", row.names = FALSE)
#304 genes left
#sort by direction - are they plus or minus LFC compared to T0?
LPos = LL_filtered %>% filter(mean_LFC > 0) %>% pull(Gene)
#14 positives - enriched relative to T0
LNeg = LL_filtered %>% filter(mean_LFC < 0) %>% pull(Gene)
#290 negatives - depleted relative to T0
#this means all of those uninformative genes missing in neurons were in the depletion list as expected
write.csv(data.frame(LPos), "FINAL_LLhitspos.csv", row.names = FALSE)
write.csv(data.frame(LNeg), "FINAL_LLhitsneg.csv", row.names = FALSE)

ULpos = read.csv("FINAL_ULhitspos.csv")
ULneg = read.csv("FINAL_ULhitsneg.csv")
LLneg = as.data.frame(LNeg)
LLpos = as.data.frame(LPos)

ss = intersect(LLneg$LNeg, ULneg$UNeg)
#17 overlapping depletions
#discard these
print(ss)
# "COX8A"  "CTH"    "DLAT"   "EHHADH" "EPRS"   "ITK"    "MMP27"  "NDUFA7" "PPP1R8" "PSAP"   "SCN2B" 
# "SI"     "TFPI"   "TUB"    "UQCRB"  "VWF"    "WARS"


#Count number of positive values (across guides) per gene (aggregated across samples)
genePos <- counts_filtered %>%
  group_by(Gene) %>%
  summarise(
    T0 = sum(T0 > 0),
    D7 = sum(D7 > 0),
    UL = sum(UL > 0),
    LL = sum(LL > 0),
    .groups = "drop"
  )
table(as.matrix(genePos[, -1]))
#   0    1    2    3    4    5    6    7    8    9   10   13   14   24   76  128  159  246 
#607 1287 1627 1672 1412 2180   94   62   71   66  170    2    1    1    1    1    1    1 

xremove = read.csv("survivals_dropped.csv")
#remove those genes entirely
counts_filtered <- counts_filtered %>%
  filter(!Gene %in% xremove$x)


#lets look at samples
gene_guides_long <- counts_filtered %>%
  pivot_longer(
    cols = c(T0, D7, UL, LL),
    names_to = "sample",
    values_to = "count"
  ) %>%
  group_by(Gene, sample) %>%
  summarise(
    nonzero_guides = sum(count > 0, na.rm = TRUE),
    total_guides = n(),
    .groups = "drop"
  )
ggplot(gene_guides_long, aes(x = sample, y = nonzero_guides)) +
  geom_boxplot() +
  labs(
    title = "Non-zero guides per gene across samples",
    x = "Sample",
    y = "Non-zero guides per gene"
  ) +
  theme_minimal()


plot_df <- gene_guides_long %>%
  mutate(category = case_when(
    nonzero_guides == 0 ~ "0",
    nonzero_guides == 1 ~ "1",
    nonzero_guides == 2 ~ "2",
    nonzero_guides == 3 ~ "3",
    nonzero_guides == 4 ~ "4",
    nonzero_guides >= 5 ~ "5+"
  ))
plot_df$category <- factor(plot_df$category, levels = c("0","1","2","3","4","5+"))

ggplot(data = plot_df, aes(x = sample, fill = category)) +
  geom_bar(position = "fill") +
  labs(
    title = "Proportion of genes by guide representation",
    y = "Proportion",
    fill = "Guide count"
  ) +
  theme_minimal()

library(scales)
library(stringr)
plot_df <- gene_guides_long %>%
  mutate(category = case_when(
    nonzero_guides == 0 ~ "0",
    nonzero_guides == 1 ~ "1",
    nonzero_guides == 2 ~ "2",
    nonzero_guides == 3 ~ "3",
    nonzero_guides == 4 ~ "4",
    nonzero_guides >= 5 ~ "5+"
  )) %>%
  count(sample, category) %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n))

plot_df$category <- factor(plot_df$category, levels = c("0","1","2","3","4","5+"))

ggplot(plot_df, aes(x = sample, y = prop, fill = category)) +
  geom_col() +
  geom_text(
    aes(label = percent(prop, accuracy = 1)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Guide representation per gene",
    y = "Percentage of genes",
    fill = "Non-zero guides"
  ) +
  theme_minimal()
#now remove non-targeting controls to make things cleaner
counts_flagged <- counts_filtered %>%
  mutate(
    neg_control = str_detect(str_to_lower(Gene), "negative_control")
  )
gene_guides_long_noNeg <- counts_flagged %>%
  filter(!neg_control) %>%
  pivot_longer(
    cols = c(T0, D7, UL, LL),
    names_to = "sample",
    values_to = "count"
  ) %>%
  group_by(Gene, sample) %>%
  summarise(
    nonzero_guides = sum(count > 0, na.rm = TRUE),
    total_guides = n(),
    .groups = "drop"
  )
plot_df_noNeg <- gene_guides_long_noNeg %>%
  mutate(category = case_when(
    nonzero_guides == 0 ~ "0",
    nonzero_guides == 1 ~ "1",
    nonzero_guides == 2 ~ "2",
    nonzero_guides == 3 ~ "3",
    nonzero_guides == 4 ~ "4",
    nonzero_guides >= 5 ~ "5+"
  )) %>%
  count(sample, category) %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n))

plot_df_noNeg$category <- factor(plot_df_noNeg$category,
                                 levels = c("0", "1", "2", "3", "4", "5+"))

ggplot(plot_df_noNeg, aes(x = sample, y = prop, fill = category)) +
  geom_col() +
  geom_text(
    aes(label = percent(prop, accuracy = 1)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Guide representation per gene (negative controls removed)",
    x = "Sample",
    y = "Percentage of genes",
    fill = "Non-zero guides"
  ) +
  theme_minimal() + coord_flip()

plot_df_noNeg <- gene_guides_long_noNeg %>%
  mutate(category = case_when(
    nonzero_guides == 0 ~ "0",
    nonzero_guides == 1 ~ "1",
    nonzero_guides == 2 ~ "2",
    nonzero_guides >= 3 ~ "3+")) %>%
  count(sample, category) %>%
  group_by(sample) %>%
  mutate(prop = n / sum(n))

plot_df_noNeg$category <- factor(plot_df_noNeg$category,
                                 levels = c("0", "1", "2", "3+"))

ggplot(plot_df_noNeg, aes(x = sample, y = prop, fill = category)) +
  geom_col() +
  geom_text(
    aes(label = percent(prop, accuracy = 1)),
    position = position_stack(vjust = 0.5),
    size = 3
  ) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    title = "Guide representation per gene (negative controls removed)",
    x = "Sample",
    y = "Percentage of genes",
    fill = "Non-zero guides"
  ) +
  theme_minimal() + coord_flip()
