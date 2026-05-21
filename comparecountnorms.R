#comparing norms
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/")
rm (list = ls())
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(readr)
#use consistent colours
sample_cols <- c(T0 = "black", D7 = "#2f95dc", UL = "#5abf41", LL = "#d95f75")
#Rug row positions
sample_levels <- c("T0", "D7", "UL", "LL")
rug_y <- c(T0 = -0.08, D7 = -0.18, UL = -0.28, LL = -0.38)
#read in data
raw_counts <- read.delim("allnonorm.count.txt", stringsAsFactors = FALSE)
ctrl_counts <- read.delim("all_counts_controlnorm/all_controlnorm.count_normalized.txt",
  stringsAsFactors = FALSE)
ctrl_counts <- ctrl_counts %>% rename(UL = "PSU", LL = "PSL")
med_counts <- read.delim("count_all/all.count_normalized.txt",
  stringsAsFactors = FALSE)
med_counts <- med_counts %>% rename(UL = "PSU", LL = "PSL")
sg_counts <- read.csv("sgcounts/allsgcounts_withsgREF.csv")
sg_counts <- sg_counts %>% rename(sgRNA = "Guide", T0 = "T0_count", D7 = "D7_count", UL = "PSU_count", LL = "PSL_count")
numeric_cols <- c("T0", "D7", "UL", "LL")

prepare_counts_for_plot <- function(df, dataset_name, remove_total_dropout = TRUE) {df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)
  if (remove_total_dropout) {
    df <- df %>%
      filter(!(T0 == 0 & D7 == 0 & UL == 0 & LL == 0))}
  
  # Replace remaining zeros with 1 before log2 transform - this proved problematic for normalised data
#instead, add pseudocount to everything
  pseudocount <- 1
    
  plot_df <- df %>%
    select(sgRNA, Gene, all_of(numeric_cols)) %>%
    pivot_longer(
      cols = all_of(numeric_cols),
      names_to = "sample",
      values_to = "count"
    ) %>%
    mutate(
      log2_count = log2(count + pseudocount),
      sample = factor(sample, levels = sample_levels),
      norm_type = dataset_name
    )
  
  neg_df <- df %>%
    filter(Gene == "negative_control") %>%
    select(sgRNA, all_of(numeric_cols)) %>%
    pivot_longer(
      cols = all_of(numeric_cols),
      names_to = "sample",
      values_to = "count"
    ) %>%
    mutate(
      log2_count = log2(count + pseudocount),
      sample = factor(sample, levels = sample_levels),
      rug_y = rug_y[as.character(sample)],
      norm_type = dataset_name
    )
  
  neg_medians <- neg_df %>%
    group_by(sample) %>%
    summarise(median_log2 = median(log2_count, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      median_log2 = round(median_log2, 2),
      norm_type = dataset_name
    )
  
  list(
    plot_df = plot_df,
    neg_df = neg_df,
    neg_medians = neg_medians
  )
}

raw_prepped  <- prepare_counts_for_plot(raw_counts,  "No normalisation", remove_total_dropout = TRUE)
ctrl_prepped <- prepare_counts_for_plot(ctrl_counts, "Negative-control normalisation", remove_total_dropout = TRUE)
med_prepped  <- prepare_counts_for_plot(med_counts,  "Median normalisation", remove_total_dropout = TRUE)
sg_prepped <- prepare_counts_for_plot(sg_counts,  "Independent counting method", remove_total_dropout = TRUE)

all_plot_df <- bind_rows(
  raw_prepped$plot_df,
  med_prepped$plot_df,
  ctrl_prepped$plot_df,
  sg_prepped$plot_df
)

all_neg_df <- bind_rows(
  raw_prepped$neg_df,
  med_prepped$neg_df,
  ctrl_prepped$neg_df,
  sg_prepped$neg_df
)

#want them on same scale/axis limits
xmin <- floor(min(all_plot_df$log2_count, na.rm = TRUE)) - 2
xmax <- ceiling(max(all_plot_df$log2_count, na.rm = TRUE)) + 2

# Plotting function
make_density_plot <- function(prepped, title_text) {
  
  plot_df <- prepped$plot_df
  neg_df <- prepped$neg_df
  
  x_right <- xmax - 1.5
  
  label_df <- data.frame(
    sample = factor(sample_levels, levels = sample_levels),
    x = x_right,
    y = rug_y[sample_levels]
  )
  
  ggplot(plot_df, aes(x = log2_count, colour = sample)) + geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.4) +
    geom_density(linewidth = 1.1, show.legend = FALSE) +
    geom_segment(
      data = neg_df,
      aes(
        x = log2_count,
        xend = log2_count,
        y = rug_y,
        yend = rug_y + 0.022, colour = sample
      ),
      inherit.aes = FALSE,
      linewidth = 0.5,
      alpha = 0.9,
      show.legend = FALSE
    ) +
    geom_text(
      data = label_df,
      aes(x = x, y = y, label = sample, colour = sample),
      inherit.aes = FALSE,
      hjust = 0,
      size = 5,
      fontface = "bold",
      show.legend = FALSE
    ) +
    scale_colour_manual(values = sample_cols, drop = FALSE) +
    scale_x_continuous(
      limits = c(xmin, xmax),
      expand = c(0, 0)) +
    scale_y_continuous(
      breaks = c(0, 0.25, 0.5),
      labels = c("0", "0.25", "0.50"),
      expand = c(0, 0)) +
    coord_cartesian(
      ylim = c(-0.5, 0.5),
      clip = "off") +
    labs(
      title = title_text,
      x = "log2(count)",
      y = "Density") +
    theme_classic(base_size = 16) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      plot.margin = margin(15, 50, 15, 15))
}
#plot
p_raw  <- make_density_plot(raw_prepped,  "No normalisation")
p_ctrl <- make_density_plot(ctrl_prepped, "Control normalisation")
p_med  <- make_density_plot(med_prepped,  "Median normalisation")
p_sgcount <- make_density_plot(sg_prepped,  "Independent count method")
# View individually
p_raw
p_ctrl
p_med
p_sgcount

bneg_df <- raw_prepped$neg_df
range(bneg_df$log2_count[bneg_df$sample == "LL"])
neg_medians <- bneg_df %>%
  group_by(sample) %>%
  summarise(median_log2 = median(log2_count), .groups = "drop")
neg_medians %>%
  mutate(median_log2 = round(median_log2, 2))
#1 T0           10.1 
# 2 D7           10.2 
# 3 UL            9.59
# 4 LL            1.58

cneg_df <- ctrl_prepped$neg_df
neg_medians <- cneg_df %>%
  group_by(sample) %>%
  summarise(median_log2 = median(log2_count), .groups = "drop")
neg_medians %>%
  mutate(median_log2 = round(median_log2, 2))
# T0            7.25
# D7            9.76
# UL           10.5 
# LL            6.72
mneg_df <- med_prepped$neg_df
neg_medians <- mneg_df %>%
  group_by(sample) %>%
  summarise(median_log2 = median(log2_count), .groups = "drop")
neg_medians %>%
  mutate(median_log2 = round(median_log2, 2))
# T0            7.19
# D7            9.62
# UL           11.8 
# LL            6.7 


#plot together:
ggplot(all_plot_df, aes(x = log2_count, colour = sample)) +
  geom_density(linewidth = 1, show.legend = TRUE) +
  geom_segment(
    data = all_neg_df,
    aes(
      x = log2_count,
      xend = log2_count,
      y = rug_y,
      yend = rug_y + 0.02, colour = sample
    ),
    inherit.aes = FALSE,
    linewidth = 0.4,
    alpha = 0.75,
    show.legend = FALSE
  ) +
  facet_wrap(~ norm_type, ncol = 2) +
  scale_colour_manual(values = sample_cols) +
  coord_cartesian(
    xlim = c(xmin, xmax),
    ylim = c(-0.5, NA),
    clip = "on"
  ) +
  labs(
    x = "log2(count)",
    y = "Density",
    colour = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
    strip.text = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 11)
  )
