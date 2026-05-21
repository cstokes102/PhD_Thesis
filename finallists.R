#making final hit list
rm (list = ls())
library(dplyr)
library(VennDiagram)
library(tidyverse)
library(ggplot2)
library(stringr)
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/sgcounts/")
PSLhits = read.csv("PSLhitconcordance.csv")
#879 hits
PSLhits$methods <- rowSums(PSLhits[,-1])
PSLhits <- PSLhits %>%
  filter(methods >= 3)
#318 if >2 and 164 if >3
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/mageck_nest/magecknest_test/")
PSUhits = read.csv("PSUmethodconco.csv")
#162
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/mageck_nest/magecknest_test/")
nesttest <- read.delim("magecknest_test_PPI_False_outliers_removal_False.gene_summary.txt", stringsAsFactors = F)
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/all_counts_controlnorm/")
ctrlnormPSU = read.csv("CS_PSUT0manLFC.csv")
psu_dir <- ctrlnormPSU %>%
  select(Gene, mean_LFC) %>%
  distinct() %>%
  mutate(PSU_direction = case_when(
    mean_LFC > 0 ~ "Positive",
    mean_LFC < 0 ~ "Negative",
    TRUE ~ "Zero"
  ))
PSUfinal <- PSUhits %>%
  left_join(psu_dir, by = "Gene") %>%
  filter(PSU_direction == "Positive")
#59
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/all_counts_controlnorm/")
ctrlnormPSL = read.csv("CS_PSLT0manLFC.csv")
psl_dir <- ctrlnormPSL %>%
  select(Gene, mean_LFC) %>%
  distinct() %>%
  mutate(PSL_direction = case_when(
    mean_LFC < 0 ~ "Negative",
    mean_LFC > 0 ~ "Positive",
    TRUE ~ "Zero"
  ))
PSLfinal <- PSLhits %>%
  left_join(psl_dir, by = "Gene") %>%
  filter(PSL_direction == "Negative")
#140
sum(is.na(PSLfinal$mean_LFC))
sum(is.na(PSUfinal$mean_LFC))
#now add in info
#NEST scores
nest_scores <- nesttest %>%
  select(Gene, Pos.beta)
#for psu
#ctrl-normalised manual LFC
ctrlnorm_scores <- ctrlnormPSU %>%
  select(Gene, mean_LFC) %>%
  rename(ctrlnorm_LFC = mean_LFC)
# non-normalised manual LFC
nonorm_scores <- read.csv("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/newCS_PSU_sigs.csv") %>%
  select(Gene, mean_LFC) %>%
  rename(nonorm_LFC = mean_LFC)
# sgcount gene scores
sg_scores <- read.csv("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/sgcounts/CS_PSU_sigs_sgcount.csv") %>%
  select(Gene, mean_LFC) %>%
  rename(sg_LFC = mean_LFC)
PSUscores <- nest_scores %>%
  left_join(ctrlnorm_scores, by="Gene") %>%
  left_join(nonorm_scores, by="Gene") %>%
  left_join(sg_scores, by="Gene")
PSUscores <- PSUfinal %>%
  left_join(PSUscores, by="Gene")
write.csv(PSUscores, "finalPSUlist.csv", row.names = F)
cor_matrix <- PSUscores %>%
  select(Pos.beta, ctrlnorm_LFC, nonorm_LFC, sg_LFC) %>%
  cor(use="pairwise.complete.obs")
cor_matrix
# Pos.beta ctrlnorm_LFC nonorm_LFC    sg_LFC
# Pos.beta     1.0000000    0.8613439  0.9688440 0.9743398
# ctrlnorm_LFC 0.8613439    1.0000000  1.0000000 0.9923107
# nonorm_LFC   0.9688440    1.0000000  1.0000000 0.9976801
# sg_LFC       0.9743398    0.9923107  0.9976801 1.0000000
#for psl
#ctrl-normalised manual LFC
ctrlnorm_scores <- ctrlnormPSL %>%
  select(Gene, mean_LFC) %>%
  rename(ctrlnorm_LFC = mean_LFC)
#non-normalised manual LFC
nonorm_scores <- read.csv("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/newCS_PSL_sigs.csv") %>%
  select(Gene, mean_LFC) %>%
  rename(nonorm_LFC = mean_LFC)
#sgcount gene scores
sg_scores <- read.csv("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/sgcounts/CS_PSL_sigs_sgcounts.csv") %>%
  select(Gene, mean_LFC) %>%
  rename(sg_LFC = mean_LFC)
nest_scores <- nesttest %>%
  select(Gene, Neg.beta)
PSLscores <- nest_scores %>%
  left_join(ctrlnorm_scores, by="Gene") %>%
  left_join(nonorm_scores, by="Gene") %>%
  left_join(sg_scores, by="Gene")
PSLscores <- PSLfinal %>%
  left_join(PSLscores, by="Gene")
#140
write.csv(PSLscores, "finalPSLlist.csv", row.names = F)
cor_matrix <- PSLscores %>%
  select(Neg.beta, ctrlnorm_LFC, nonorm_LFC, sg_LFC) %>%
  cor(use="pairwise.complete.obs")
cor_matrix
# Neg.beta ctrlnorm_LFC nonorm_LFC    sg_LFC
# Neg.beta     1.0000000    0.1332053  0.1332053 0.1545303
# ctrlnorm_LFC 0.1332053    1.0000000  1.0000000 0.9581731
# nonorm_LFC   0.1332053    1.0000000  1.0000000 0.9581731
# sg_LFC       0.1545303    0.9581731  0.9581731 1.0000000
#nest is on a completely different scale

