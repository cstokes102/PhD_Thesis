#method comparison scatter plots
rm (list = ls())
library(dplyr)
library(VennDiagram)
library(tidyverse)
library(ggplot2)
library(stringr)
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/")
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
sg = read.csv("sgcounts/FINAL_ULhits.csv")
sg <- anti_join(sg, genes_toremove, by = "Gene")
#133
sgUPos <- sg %>% filter(mean_LFC > 0) %>% pull(Gene)
#61
sgUNeg <- sg %>% filter(mean_LFC < 0) %>% pull(Gene)
#72
nonorm = read.csv("FINAL_allnonormhits.csv")
nonorm = nonorm %>%
  filter(abs(nonorm$Zscore) > 1.5)
nUPos <- nonorm %>% filter(mean_LFC > 0) %>% pull(Gene)
#66
nUNeg <- nonorm %>% filter(mean_LFC < 0) %>% pull(Gene)
#65 genes
ctrlnorm = read.csv("all_counts_controlnorm/FINAL_ULhits.csv")
#142 genes
ctrlnorm <- anti_join(ctrlnorm, genes_toremove, by = "Gene")
#131
cUPos <- ctrlnorm %>% filter(mean_LFC > 0) %>% pull(Gene)
#71
cUNeg <- ctrlnorm %>% filter(mean_LFC < 0) %>% pull(Gene)
#60
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/mageck_nest/")
PSUT0nest <- read.delim("nest_PSUT0/nest_PSUT0_PPI_False_outliers_removal_False.gene_summary.txt", stringsAsFactors = F)
#nest -k all.count.txt -d pairmatrix.txt -n nest_PSUT0 --norm-method control -e negative_control 
validPSUT0 <- PSUT0nest %>% 
  filter(!is.na(Pos.z) & is.finite(Pos.z))
#remove genes with fewer than 3 guides
validPSUT0 <- validPSUT0 %>%
  filter(sgRNA >= 3) %>% rename(beta_PSUT0 = Pos.beta)
#2288 genes --> 2288 genes
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/mageck_nest/magecknest_test/")
nesttest <- read.delim("magecknest_test_PPI_False_outliers_removal_False.gene_summary.txt", stringsAsFactors = F)
#nest -k ../all_counts_controlnorm/all_controlnorm.count.txt -d new_designmatrix.txt --norm-method control -e negative_control 
validtest <- nesttest %>% 
  filter(!is.na(Pos.z) & is.finite(Pos.z))
#has 1865 genes in it
#keep only useful columns from NEST
validtest <- nesttest %>% 
  filter(!is.na(Pos.z) & is.finite(Pos.z)) %>%
  select(Gene, sgRNA,
         Pos.beta, Pos.z, Pos.wald.p.value, Pos.wald.fdr)
#remove genes with fewer than 3 guides
validtest <- validtest %>%
  filter(sgRNA >= 3) %>% rename(beta_test = Pos.beta)
#1849 genes now
validtest <- anti_join(validtest, genes_toremove, by = "Gene")
#1710
validPSUT0 <- anti_join(validPSUT0, genes_toremove, by = "Gene")
#2145
#merge betas
merged_beta <- validtest %>%
  select(Gene, beta_test) %>%
  inner_join(validPSUT0 %>% select(Gene, beta_PSUT0), by = "Gene")

merged_beta <- merged_beta %>%
  mutate(
    hit_direction = case_when(
      Gene %in% nUPos ~ "Positive",
      Gene %in% nUNeg ~ "Negative",
      TRUE ~ "Non-hit"
    )
  )

merged_beta$hit_direction <- factor(
  merged_beta$hit_direction,
  levels = c("Non-hit", "Negative", "Positive"))
#add in nonorm hits
combined <- merged_beta %>%
  left_join(nonorm, by = "Gene")
hit_cols <- c(
  "Non-hit" = "grey70",
  "Negative" = "#f26f63",
  "Positive" = "#15b8c4")
base_theme <- theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13))
#compare nests
p1 <- ggplot(merged_beta, aes(x = beta_PSUT0, y = beta_test)) +
  geom_point(
    data = subset(merged_beta, hit_direction == "Non-hit"),
    colour = "grey70",
    alpha = 0.35,
    size = 2.5
  ) +
  geom_point(
    data = subset(merged_beta, hit_direction != "Non-hit"),
    aes(colour = hit_direction),
    size = 4
  ) +
  scale_colour_manual(
    values = hit_cols[c("Negative", "Positive")],
    breaks = c("Negative", "Positive")
  ) +
  labs(
    x = "beta_PSUT0",
    y = "beta_test",
    colour = "Hit direction"
  ) + base_theme
p1
p2 <- ggplot(combined, aes(x = mean_LFC, y = beta_test)) +
  geom_point(
    data = subset(combined, hit_direction == "Non-hit"),
    colour = "black",
    alpha = 0.35,
    size = 2.5
  ) +
  geom_point(
    data = subset(combined, hit_direction != "Non-hit"),
    aes(colour = hit_direction),
    size = 4
  ) +
  scale_colour_manual(
    values = hit_cols[c("Negative", "Positive")],
    breaks = c("Negative", "Positive")
  ) +
  labs(
    x = "Manual mean LFC",
    y = "MAGeCK-NEST beta",
    colour = "PSU_direction"
  ) +
  base_theme
p2

UL_manual_scores <- nonorm %>%
  select(Gene, mean_LFC, sd_LFC, n_guides, Zscore) %>%
  rename(manual_LFC = mean_LFC,
         manual_sd = sd_LFC,
         manual_n_guides = n_guides,
         manual_Z = Zscore)
UL_pos_hits <- data.frame(Gene = nUPos)
#manual LFC v beta scatterplot
UL_compare <- UL_manual_scores %>%
  left_join(validtest %>% select(Gene, beta_test, Pos.wald.fdr), by = "Gene")
p_UL_scatter <- ggplot(UL_compare, aes(manual_LFC, beta_test)) +
  geom_point(alpha = 0.15, color = "grey50") +
  geom_point(data = UL_compare %>% filter(Gene %in% nUPos),
             size = 2.5) +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(
    title = "UL: manual LFC vs MAGeCK-NEST beta",
    x = "Manual mean LFC (UL vs T0)",
    y = "NEST Pos.beta"
  )
p_UL_scatter

UL_rank <- UL_compare %>%
  filter(Gene %in% nUPos, !is.na(beta_test), !is.na(manual_LFC)) %>%
  arrange(desc(beta_test)) %>%
  mutate(nest_rank = row_number()) %>%
  arrange(desc(manual_LFC)) %>%
  mutate(manual_rank = row_number())
cor(UL_rank$manual_rank, UL_rank$nest_rank,
    method = "spearman", use = "complete.obs")
#0.7912754
LL_rank <- LL_compare %>%
  filter(Gene %in% LNeg, !is.na(Neg.beta), !is.na(manual_LFC)) %>%
  arrange(Neg.beta) %>%
  mutate(nest_rank = row_number()) %>%
  arrange(manual_LFC) %>%
  mutate(manual_rank = row_number())
cor(LL_rank$manual_rank, LL_rank$nest_rank,
    method = "spearman", use = "complete.obs")
#0.5422337
p_UL_rank <- ggplot(UL_rank, aes(manual_rank, nest_rank)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(
    title = "UL hit rank concordance",
    x = "Manual rank",
    y = "NEST rank"
  )
p_UL_rank
p_LL_rank <- ggplot(LL_rank, aes(manual_rank, nest_rank)) +
  geom_point() +
  geom_smooth(method = "lm") +
  theme_classic() +
  labs(
    title = "LL hit rank concordance",
    x = "Manual rank",
    y = "NEST rank"
  )

p_UL_rank
p_LL_rank
# UL method scores
UL_sg_scores <- read.csv("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/sgcounts/FINAL_ULhits.csv") %>%
  select(Gene, mean_LFC) %>%
  rename(sgcount = mean_LFC)
UL_ctrlnorm_scores <- read.csv("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/all_counts_controlnorm/FINAL_ULhits.csv") %>%
  select(Gene, mean_LFC) %>%
  rename(ctrlnorm = mean_LFC)
# LL method scores
LL_sg_scores <- read.csv("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/sgcounts/FINAL_LLhits.csv") %>%
  select(Gene, mean_LFC) %>%
  rename(sgcount = mean_LFC)
LL_ctrlnorm_scores <- read.csv("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/all_counts_controlnorm/FINAL_LLhits.csv") %>%
  select(Gene, mean_LFC) %>%
  rename(ctrlnorm = mean_LFC)
UL_heat <- UL_pos_hits %>%
  left_join(nest_valid %>% select(Gene, Pos.beta) %>% rename(NEST = Pos.beta), by = "Gene") %>%
  left_join(UL_manual_scores %>% select(Gene, manual_LFC) %>% rename(manual = manual_LFC), by = "Gene") %>%
  left_join(UL_sg_scores, by = "Gene") %>%
  left_join(UL_ctrlnorm_scores, by = "Gene")
UL_mat <- UL_heat %>%
  select(Gene, NEST, manual, sgcount, ctrlnorm) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  drop_na() %>%
  column_to_rownames("Gene") %>%
  as.matrix()
UL_mat_scaled <- scale(UL_mat)
pheatmap(
  UL_mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Positive UL hits: concordance across methods",
  fontsize_row = 6
)
UL_order <- UL_compare %>%
  filter(Gene %in% UPos, !is.na(Pos.beta)) %>%
  arrange(desc(Pos.beta)) %>%
  select(Gene)
UL_heat2 <- UL_order %>%
  left_join(nest_valid %>% select(Gene, Pos.beta) %>% rename(NEST = Pos.beta), by = "Gene") %>%
  left_join(UL_manual_scores %>% select(Gene, manual_LFC) %>% rename(manual = manual_LFC), by = "Gene") %>%
  left_join(UL_sg_scores, by = "Gene") %>%
  left_join(UL_ctrlnorm_scores, by = "Gene")
UL_mat2 <- UL_heat2 %>%
  select(Gene, NEST, manual, sgcount, ctrlnorm) %>%
  drop_na() %>%
  column_to_rownames("Gene") %>%
  as.matrix()
UL_mat2_scaled <- scale(UL_mat2)
pheatmap(
  UL_mat2_scaled,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Positive UL hits ordered by NEST rank",
  fontsize_row = 6
)
cor(UL_mat2, method = "spearman", use = "pairwise.complete.obs")
#NEST    manual   sgcount  ctrlnorm
#NEST     1.0000000 0.8218823 0.8324475 0.8253176
#manual   0.8218823 1.0000000 0.9853513 0.9994815
#sgcount  0.8324475 0.9853513 1.0000000 0.9857402
#ctrlnorm 0.8253176 0.9994815 0.9857402 1.0000000



LL_heat <- LL_neg_hits %>%
  left_join(nest_valid %>% select(Gene, Neg.beta) %>% rename(NEST = Neg.beta), by = "Gene") %>%
  left_join(LL_manual_scores %>% select(Gene, manual_LFC) %>% rename(manual = manual_LFC), by = "Gene") %>%
  left_join(LL_sg_scores, by = "Gene") %>%
  left_join(LL_ctrlnorm_scores, by = "Gene")

LL_mat <- LL_heat %>%
  select(Gene, NEST, manual, sgcount, ctrlnorm) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  drop_na() %>%
  column_to_rownames("Gene") %>%
  as.matrix()

LL_mat_scaled <- scale(LL_mat)

pheatmap(
  LL_mat_scaled,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  main = "Negative LL hits: concordance across methods",
  fontsize_row = 6
)


