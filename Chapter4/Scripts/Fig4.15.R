#nonorm guide representation
rm (list = ls())
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(scales)
library(stringr)
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/")
#read in data
raw_counts <- read.delim("allnonorm.count.txt", stringsAsFactors = FALSE)
sample_levels <- c("T0", "D7", "UL", "LL")
numeric_cols <- c("T0", "D7", "UL", "LL")
raw_counts[numeric_cols] <- lapply(raw_counts[numeric_cols], function(x) as.numeric(as.character(x)))
#how many guides exist for each gene in the library
guides_per_gene <- raw_counts %>%
  group_by(Gene) %>%
  summarise(total_guides = n(), .groups = "drop")
#look for total guide dropout
totaldropout <- raw_counts %>%
  filter(T0 == 0 & D7 == 0 & UL == 0 & LL == 0)
#433 dropouts
counts_filtered <- anti_join(raw_counts, totaldropout, by = "sgRNA")
#12575 guides
#lets look at samples
gene_guides_long <- counts_filtered %>%
  pivot_longer(
    cols = all_of(sample_levels),
    names_to = "sample",
    values_to = "count") %>%
  group_by(Gene, sample) %>%
  summarise(
    nonzero_guides = sum(count > 0, na.rm = TRUE),
    total_guides = n(),
    .groups = "drop") %>%
  mutate(sample = factor(sample, levels = sample_levels))
ggplot(gene_guides_long, aes(x = sample, y = nonzero_guides)) +
  geom_boxplot() +
  labs(
    title = "Non-zero guides per gene across samples",
    x = "Sample",
    y = "Non-zero guides per gene") +
  theme_minimal()
#plot guides per gene
plot_df <- gene_guides_long %>%
  mutate(category = case_when(
    nonzero_guides == 0 ~ "0",
    nonzero_guides == 1 ~ "1",
    nonzero_guides == 2 ~ "2",
    nonzero_guides >= 3 ~ "3+" ))
plot_df$category <- factor(plot_df$category, levels = c("0","1","2","3+"))
ggplot(data = plot_df, aes(x = sample, fill = category)) +
  geom_bar(position = "fill") +
  labs(
    title = "Non-zero guides per gene",
    y = "Proportion",
    fill = "Guide count") +
  theme_minimal()

#now treat guides with <10 reads as absent
raw_counts <- raw_counts %>%
  mutate(across(all_of(numeric_cols), ~ ifelse(. < 10, 0, .)))
#how many guides exist for each gene in the library
guides_per_gene <- raw_counts %>%
  group_by(Gene) %>%
  summarise(total_guides = n(), .groups = "drop")
#remove guides absent in all samples
totaldropout <- raw_counts %>%
  filter(T0 == 0 & D7 == 0 & UL == 0 & LL == 0)
counts_filtered <- anti_join(raw_counts, totaldropout, by = "sgRNA")
# summarise effective guide representation per gene per sample
gene_guides_long <- counts_filtered %>%
  pivot_longer(
    cols = all_of(sample_levels),
    names_to = "sample",
    values_to = "count"
  ) %>%
  group_by(Gene, sample) %>%
  summarise(
    nonzero_guides = sum(count > 0, na.rm = TRUE),
    total_guides = n(),
    .groups = "drop"
  ) %>%
  mutate(sample = factor(sample, levels = sample_levels))
# categorise genes by number of sufficiently represented guides
plot_df <- gene_guides_long %>%
  mutate(
    category = case_when(
      nonzero_guides == 0 ~ "0",
      nonzero_guides == 1 ~ "1",
      nonzero_guides == 2 ~ "2",
      nonzero_guides >= 3 ~ "3+"
    ),
    category = factor(category, levels = c("0", "1", "2", "3+"))
  ) %>%
  count(sample, category, name = "n_genes") %>%
  group_by(sample) %>%
  mutate(percent_genes = 100 * n_genes / sum(n_genes)) %>%
  ungroup()
# stacked barplot
ggplot(plot_df, aes(x = sample, y = percent_genes, fill = category)) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = ifelse(percent_genes >= 4, paste0(round(percent_genes), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    colour = "black"
  ) +
  scale_fill_manual(
    values = c(
      "0"  = "#f08080",   # dropout
      "1"  = "#fcb85e",   # weak
      "2"  = "#fdd55a",   # moderate
      "3+" = "#95c11f"    # strong
    ),
    name = "Guides per gene\nwith \u226510 reads"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    title = "Guide representation per gene",
    x = "Sample",
    y = "% of library genes"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )


#This now shows effective representation, not just non-zero representation.
#% of library genes here means the percentage of genes in the filtered library, after removing guides absent in all samples.
#Guide representation per gene after applying a minimum 10-read threshold
#now remove non-targeting controls to make things cleaner
counts_flagged <- counts_filtered %>%
  mutate(
    neg_control = str_detect(str_to_lower(Gene), "negative_control"))
gene_guides_long_noNeg <- counts_flagged %>% filter(!neg_control) %>%
  pivot_longer(
    cols = all_of(sample_levels),
    names_to = "sample",
    values_to = "count"
  ) %>%
  group_by(Gene, sample) %>%
  summarise(
    nonzero_guides = sum(count > 0, na.rm = TRUE),
    total_guides = n(),
    .groups = "drop"
  ) %>%
  mutate(sample = factor(sample, levels = sample_levels))
# categorise genes by number of sufficiently represented guides
plot_df1 <- gene_guides_long_noNeg %>%
  mutate(
    category = case_when(
      nonzero_guides == 0 ~ "0",
      nonzero_guides == 1 ~ "1",
      nonzero_guides == 2 ~ "2",
      nonzero_guides >= 3 ~ "3+"
    ),
    category = factor(category, levels = c("0", "1", "2", "3+"))
  ) %>%
  count(sample, category, name = "n_genes") %>%
  group_by(sample) %>%
  mutate(percent_genes = 100 * n_genes / sum(n_genes)) %>%
  ungroup()
# stacked barplot
ggplot(plot_df1, aes(x = sample, y = percent_genes, fill = category)) +
  geom_col(width = 0.75) +
  geom_text(
    aes(label = ifelse(percent_genes >= 4, paste0(round(percent_genes), "%"), "")),
    position = position_stack(vjust = 0.5),
    size = 3,
    colour = "black"
  ) +
  scale_fill_manual(
    values = c(
      "0"  = "#f08080",   # dropout
      "1"  = "#fcb85e",   # weak
      "2"  = "#fdd55a",   # moderate
      "3+" = "#95c11f"    # strong
    ),
    name = "Guides per gene\nwith \u226510 reads"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    title = "Guide representation per gene",
    x = "Sample",
    y = "% of library genes"
  ) +
  theme_classic(base_size = 15) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )
