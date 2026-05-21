dev.off()
rm (list = ls())
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/mageck_nest/magecknest_test/")
library("clusterProfiler", lib = "~/Rlibs")
library("dplyr")
library(readr)
library(ggplot2)
library(tidyr)
library(knitr)
library(limma)
library(edgeR)
nesttest <- read.delim("magecknest_test_PPI_False_outliers_removal_False.gene_summary.txt", stringsAsFactors = F)
testguides <- read.delim("magecknest_test_PPI_False_outliers_removal_False.sgrna_summary.txt", stringsAsFactors = F, header = F)
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/mageck_nest/nest_controlnormrawcounts/")
newnest  <- read.delim("nest_controlnormrawcounts_PPI_False_outliers_removal_False.gene_summary.txt", stringsAsFactors = F)
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/mageck_nest/")
nesttest <- nesttest %>%
  select(Gene, Pos.beta, Pos.z, Pos.wald.p.value, Pos.wald.fdr)
newnest <- newnest %>%
  select(Gene, Pos.beta, Pos.z, Pos.wald.p.value, Pos.wald.fdr)
#how many are NA
table(is.nan(nesttest$Pos.beta)) #false so everything does technically have a beta score
table(is.na(nesttest$Pos.z)) #false 1865 true   444 
table(is.nan(newnest$Pos.beta)) #false so everything does technically have a beta score
table(is.na(newnest$Pos.z)) #false 466 true 1843 
valid <- newnest %>% 
  filter(!is.na(Pos.z) & is.finite(Pos.z))
validtest <- nesttest %>% 
  filter(!is.na(Pos.z) & is.finite(Pos.z))
merged <- inner_join(z_df, validtest, by = "Gene")
cor.test(merged$Zscore, merged$Pos.beta, method = "spearman")  # rank correlation
# rho = 0.8328612 
#S = 174666112, p-value < 2.2e-16
cor.test(merged$Zscore, merged$Pos.beta, method = "pearson")   # linear correlation
#t = 27.275, df = 1842, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval: 0.5030317 0.5681018
#sample estimates:cor 0.5363633 
ggplot(merged, aes(x = Zscore, y = Pos.beta)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Comparison of Manual Gene-level Z-scores vs MAGeCK-NEST β-scores",
    subtitle = "PSU vs T0 (positive selection)",
    x = "Manual per-gene Z-score",
    y = "MAGeCK-NEST (test) Pos.beta")

merged1 <- inner_join(z_df, valid, by = "Gene")
cor.test(merged1$Zscore, merged1$Pos.beta, method = "spearman")  # rank correlation
#S = 3233859, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#rho =  0.7842054 
cor.test(merged1$Zscore, merged1$Pos.beta, method = "pearson")   # linear correlation
#t = 19.282, df = 446, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval: 0.6203711 0.7218181
#sample estimates:cor 0.6742632 
ggplot(merged1, aes(x = Zscore, y = Pos.beta)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Comparison of Manual Gene-level Z-scores vs MAGeCK-NEST β-scores",
    subtitle = "PSU vs T0 (positive selection)",
    x = "Manual per-gene Z-score",
    y = "MAGeCK-NEST Pos.beta")
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/mageck_nest/nest_PSUT0/")
newnest  <- read.delim("nest_PSUT0_PPI_False_outliers_removal_False.gene_summary.txt", stringsAsFactors = F)
table(is.nan(newnest$Pos.beta)) #false so everything does technically have a beta score
table(is.na(newnest$Pos.z)) #false 466 true 1843 
valid <- newnest %>% 
  filter(!is.na(Pos.z) & is.finite(Pos.z))
merged2 <- inner_join(z_df, valid, by = "Gene")
cor.test(merged2$Zscore, merged2$Pos.beta, method = "spearman")  # rank correlation
#S = 75716847, p-value < 2.2e-16
#alternative hypothesis: true rho is not equal to 0
#rho =  0.9618209  - almost perfect rank correlation!!
cor.test(merged2$Zscore, merged2$Pos.beta, method = "pearson")   # linear correlation
#t = 31.661, df = 2281, p-value < 2.2e-16
#95 percent confidence interval: 0.5233777 0.5804057
#sample estimates:cor 0.552538
#likely lower cos beta scores and z scores have diff scales (not linear relationship)
ggplot(merged2, aes(x = Zscore, y = Pos.beta)) +
  geom_point(alpha = 0.6, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Comparison of Manual Gene-level Z-scores vs MAGeCK-NEST β-scores",
    subtitle = "PSU vs T0 (positive selection)",
    x = "Manual per-gene Z-score",
    y = "MAGeCK-NEST Pos.beta")
