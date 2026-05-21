rm (list = ls())
library(stringr)
library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/mageck_nest/")
nest <- read.delim("realdeal/realdeal_PPI_False_outliers_removal_False.gene_summary.txt", stringsAsFactors = FALSE)
nest_valid <- nest %>%
  filter(!is.na(Pos.z) & is.finite(Pos.z),
         !is.na(Neg.z) & is.finite(Neg.z),!is.na(CTRL.z) & is.finite(CTRL.z),
         sgRNA >= 3)
#1094 genes
nest_contrasts <- nest_valid %>%
  filter(!is.na(Pos.beta), !is.na(CTRL.beta), sgRNA >= 3) %>%
  select(Gene, CTRL.beta, Pos.beta, Neg.beta)
nest_contrasts <- nest_contrasts %>%
  mutate(
    UL_vs_D7 = Pos.beta - CTRL.beta,
    UL_vs_mean_other = Pos.beta - ((CTRL.beta + Neg.beta) / 2),
    UL_vs_max_other = Pos.beta - pmax(CTRL.beta, Neg.beta))
# get pathways
msig <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)

#ranked vector: UL v mean of others
ranksULvmean <- nest_contrasts$UL_vs_mean_other
names(ranksULvmean) <- nest_contrasts$Gene
ranksULvmean <- sort(ranksULvmean, decreasing = TRUE)
summary(ranksULvmean)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-7.29106 -1.01031  0.01677  0.16029  1.59199  6.61715 
#show dist
hist(ranksULvmean, breaks = 50)
#gsea expects roughly symmetric rank distribution 
fgsea_vmean <- fgsea(pathways = pathways, stats = ranksULvmean, minSize = 10, maxSize = 500, nPermSimple = 2000)
# There are ties in the preranked stats (0.06% of the list).
# The order of those tied genes will be arbitrary, which may produce unexpected results.
fgsea_vmean %>%
  arrange(padj) %>%
  head(20)
#pathway         pval        padj    log2err         ES        NES  size
# HALLMARK_KRAS_SIGNALING_UP 0.0001580088 0.005846325 0.51884808 -0.7119108 -2.1194883    18
#       HALLMARK_KRAS_SIGNALING_DN 0.0324232082 0.599829352 0.23336551  0.5020850  1.5058465    22
#3:            HALLMARK_ADIPOGENESIS 0.2278325123 0.669391584 0.10119625 -0.3372020 -1.1630670    33
#4:     HALLMARK_ALLOGRAFT_REJECTION 0.1400679117 0.669391584 0.10776008  0.4189273  1.2757285    24
fgsea_vmean$leadingEdge[[1]]
plot_vm <- fgsea_vmean %>% filter(abs(NES) > 1.3) %>%
  mutate(
    pathway_clean = str_remove(pathway, "^HALLMARK_"),
    pathway_clean = str_replace_all(pathway_clean, "_", " "),
    sig_label = ifelse(padj < 0.25, "*", "")
  ) %>%
  arrange(NES)
ggplot(plot_vm, aes(x = reorder(pathway_clean, NES), y = NES)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Normalized enrichment score (NES)",
    title = "Hallmark pathway enrichment in UL (v mean of neg/d7)"
  )

ranks <- nest_contrasts$UL_vs_D7
names(ranks) <- nest_contrasts$Gene
ranks <- sort(ranks, decreasing = TRUE)
fgsea_hallmark <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 500, nPermSimple = 2000)
#There are ties in the preranked stats (0.09% of the list).
fgsea_hallmark %>%
  arrange(padj) %>%
  head(20)
fgsea_hallmark$leadingEdge[[1]]
plot_df <- fgsea_hallmark %>% filter(abs(NES) > 1.3) %>%
  mutate(
    pathway_clean = str_remove(pathway, "^HALLMARK_"),
    pathway_clean = str_replace_all(pathway_clean, "_", " "),
    sig_label = ifelse(padj < 0.25, "*", ""),
  ) %>%
  arrange(NES)
ggplot(plot_df, aes(x = reorder(pathway_clean, NES), y = NES)) +
  geom_col() +
  geom_text(aes(label = sig_label), hjust = -0.2, size = 6) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Normalized enrichment score (NES)",
    title = "Hallmark pathway enrichment UL v D7"
  )


#ranked vector: UL only
ranksUL <- nest_contrasts$Pos.beta
names(ranksUL) <- nest_contrasts$Gene
ranksUL <- sort(ranksUL, decreasing = TRUE)
summary(ranksUL)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-6.9871 -2.5779 -0.7992 -1.5018 -0.1021  6.0663 
#show dist
hist(ranksUL, breaks = 50)
fgsea_UL <- fgsea(pathways = pathways, stats = ranksUL, minSize = 10, maxSize = 500, nPermSimple = 2000)
# There are ties in the preranked stats (0.27% of the list).
# The order of those tied genes will be arbitrary, which may produce unexpected results.
fgsea_UL %>%
  arrange(padj) %>%
  head(20)
#                                     pathway        pval       padj    log2err         ES        NES  size
#  1:                  HALLMARK_G2M_CHECKPOINT 0.003384269 0.06260898 0.43170770 -0.7739088 -1.4942756    12
#2:       HALLMARK_OXIDATIVE_PHOSPHORYLATION 0.003200340 0.06260898 0.43170770 -0.6315870 -1.3397038    40
#3:         HALLMARK_CHOLESTEROL_HOMEOSTASIS 0.005846759 0.07211003 0.40701792 -0.7800350 -1.4647925    10
#4:                       HALLMARK_APOPTOSIS 0.017623364 0.13041289 0.24348488 -0.6609509 -1.3414582    20
#5:               HALLMARK_KRAS_SIGNALING_UP 0.014653866 0.13041289 0.26829569 -0.6796784 -1.3682737    18
#6:                     HALLMARK_E2F_TARGETS 0.035007610 0.21588026 0.17130086 -0.6660412 -1.3284553    16
#7:           HALLMARK_FATTY_ACID_METABOLISM 0.073536768 0.30272342 0.11474675 -0.5675672 -1.1957687    36
#12:                HALLMARK_MTORC1_SIGNALING 0.242957746 0.71300858 0.05721269 -0.5504502 -1.1296716    23
#13: HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY 0.250516529 0.71300858 0.05705686 -0.5930322 -1.1450363    12
#17:                      HALLMARK_PEROXISOME 0.439752832 0.90393638 0.03739256 -0.5371742 -1.0454315    13
#18:           HALLMARK_XENOBIOTIC_METABOLISM 0.434782609 0.90393638 0.03678762 -0.4861504 -1.0312080    40
fgsea_UL$leadingEdge[[1]]
# "C3"       "COX7B"    "CYC1"     "DLAT"     "DLD"      "SOD1"     "ACOX1"    "SLC25A10" "MAP4K3"  
# "LTC4S"    "UCK1"     "PHYH"     "MCCC1" 
plot_df <- fgsea_UL %>% filter(abs(NES) > 1.3) %>%
  mutate(
    pathway_clean = str_remove(pathway, "^HALLMARK_"),
    pathway_clean = str_replace_all(pathway_clean, "_", " "),
    sig_label = ifelse(padj < 0.25, "*", ""),
  ) %>%
  arrange(NES)
ggplot(plot_df, aes(x = reorder(pathway_clean, NES), y = NES)) +
  geom_col() +
  geom_text(aes(label = sig_label), hjust = -0.2, size = 6) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Normalized enrichment score (NES)",
    title = "Hallmark pathway enrichment in UL"
  )
