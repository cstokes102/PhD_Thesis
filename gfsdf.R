library(stringr)
genes_toremove <- dropout_vs_total %>%
  filter(guides_remaining <= 2)
#143 genes with 2 or fewer guides remaining in neurons (not informative/reliable)
#remove these genes from the hit list
#z_df has 2294 genes
UL <- z_df %>% filter(!Gene %in% genes_toremove$Gene)
#we're left with 2151 genes
write.csv(UL, "althehits.csv", row.names = FALSE)
library(fgsea)
library(msigdbr)

# get pathways
msig <- msigdbr(species = "Homo sapiens", category = "C2")
pathways <- split(msig$gene_symbol, msig$gs_name)
UL_manual_scores <- UL %>%
  select(Gene, mean_LFC, sd_LFC, n_guides, Zscore) %>%
  rename(manual_LFC = mean_LFC,
         manual_sd = sd_LFC,
         manual_n_guides = n_guides,
         manual_Z = Zscore)
# ranked vector
ranks <- UL_manual_scores$manual_LFC
names(ranks) <- UL_manual_scores$Gene
ranks <- sort(ranks, decreasing = TRUE)

fgsea_res <- fgsea(pathways = pathways, stats = ranks)
#Warning messages:
# 1: In fgseaMultilevel(pathways = pathways, stats = stats, minSize = minSize,  :
#There were 10 pathways for which P-values were not calculated properly due to unbalanced (positive and negative) gene-level statistic values. For such pathways pval, padj, NES, log2err are set to NA. You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = 10000)
summary(ranks)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-10.414  -5.732  -4.073  -4.231  -2.796   2.198
hist(ranks, breaks = 50)
#this shows me that the distribution is shifted below zero/dominated by depletion signals
#gsea expects roughly symmetric rank distribution 
#2: In fgseaMultilevel(pathways = pathways, stats = stats, minSize = minSize,  :
#For some of the pathways the P-values were likely overestimated. For such pathways log2err is set to NA.
#gonna median-centre the ranks
ranks_centered <- ranks - median(ranks)
summary(ranks_centered)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#-6.3413 -1.6592  0.0000 -0.1584  1.2770  6.2710
# look at top pathways
fgsea_res <- fgsea(pathways = pathways, stats = ranks_centered, minSize = 10, maxSize = 500)
fgsea_res %>%
  arrange(padj) %>%
  head(20)
#narrow down pathways
msig <- msigdbr(species = "Homo sapiens", category = "C2", subcollection = "CP:REACTOME")
pathways <- split(msig$gene_symbol, msig$gs_name)
fgsea_res <- fgsea(pathways = pathways, stats = ranks_centered, minSize = 10, maxSize = 500)
fgsea_res %>%
  arrange(padj) %>%
  head(20)
fgsea_res$leadingEdge[[1]]
#try hallmark pathways
msig <- msigdbr(species = "Homo sapiens", category = "H")
pathways <- split(msig$gene_symbol, msig$gs_name)
fgsea_hallmark <- fgsea(pathways = pathways, stats = ranks_centered, minSize = 10, maxSize = 500)
fgsea_hallmark %>%
  arrange(padj) %>%
  head(20)
fgsea_hallmark$leadingEdge[[1]]
#"COX8A" "C3"    "CYC1"  "ACOX1" "CPT2"  "DLAT"
collapsed <- collapsePathways(
  fgsea_res[order(pval)][1:20],
  pathways,
  stats = ranks_centered
)

mainPathways <- fgsea_res[pathway %in% collapsed$mainPathways]
msig_reactome <- msigdbr(
  species = "Homo sapiens",
  category = "C2",
  subcollection = "CP:REACTOME"
)

reactome_pathways <- split(msig_reactome$gene_symbol, msig_reactome$gs_name)

fgsea_reactome <- fgsea(
  pathways = reactome_pathways,
  stats = ranks_centered,
  minSize = 10,
  maxSize = 500,
  nPermSimple = 10000
)

reactome_subset <- fgsea_reactome %>%
  filter(str_detect(
    pathway,
    regex("AUTOPH|MITO|CHOLESTEROL|STEROL|LIPO|LYSO|LYSOSOME|GLYC", ignore_case = TRUE)
  )) %>%
  arrange(pval)

plot_df <- fgsea_hallmark %>%
  mutate(
    pathway_clean = str_remove(pathway, "^HALLMARK_"),
    pathway_clean = str_replace_all(pathway_clean, "_", " ")
  ) %>%
  arrange(NES)
ggplot(plot_df, aes(x = reorder(pathway_clean, NES), y = NES)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Normalized enrichment score (NES)",
    title = "Hallmark pathway enrichment in UL"
  )

hallmark_plot <- plot_df %>%
  mutate(direction = ifelse(NES > 0, "Positive", "Negative"))
hallmark_plot <- hallmark_plot %>%
  mutate(sig = padj < 0.25)
ggplot(hallmark_plot, aes(x = reorder(pathway_clean, NES), y = NES, fill = direction)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "NES",
    title = "Hallmark pathway enrichment in UL"
  )
hallmark_top <- hallmark_plot %>%
  slice_max(order_by = abs(NES), n = 12)

ggplot(hallmark_top, aes(x = reorder(pathway_clean, NES), y = NES, fill = direction)) +
  geom_col() +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "NES",
    title = "Top Hallmark pathway shifts in UL"
  )
hallmark_top <- hallmark_top %>%
  mutate(sig = padj < 0.25)
ggplot(hallmark_top, aes(x = reorder(pathway_clean, NES), y = NES, fill = direction)) + geom_text(aes(label = ifelse(sig, "*", "")), hjust = -0.2) +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "NES",
    title = "Top Hallmark pathways in UL"
  )
leading_df <- hallmark_top %>%
  select(pathway, leadingEdge) %>%
  tidyr::unnest(leadingEdge) %>%
  rename(Gene = leadingEdge)
gene_overlap <- leading_df %>%
  count(Gene, sort = TRUE)

#look at enriched payhways only
hallmark_pos <- hallmark_plot %>%
  filter(NES > 0)
ggplot(hallmark_pos, aes(x = reorder(pathway_clean, NES), y = NES)) +
  geom_col(fill = "#1bb3b1") +
  coord_flip() +
  theme_classic() +
  labs(
    x = NULL,
    y = "Normalized enrichment score (NES)",
    title = "Enriched Hallmark pathways in UL"
  )



ul_hits <- read.csv("FINAL_ULhitspos.csv")
ul_hits <- ul_hits$UPos
background <- UL_manual_scores$Gene
library("clusterProfiler", lib = "~/Rlibs")
library(msigdbr)

msig_h <- msigdbr(species = "Homo sapiens", category = "H")

hallmark_list <- split(msig_h$gene_symbol, msig_h$gs_name)

ora_res <- enricher(
  gene = ul_hits,
  universe = background,
  TERM2GENE = msig_h[, c("gs_name", "gene_symbol")], minGSSize = 5
)
#to make more UL-specific, I could try UL LFC - [(D7LFC + LL_LFC)/2]
