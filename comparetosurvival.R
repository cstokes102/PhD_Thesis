setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/")
rawcounts = read.delim("all_nonorm.count.txt")
#Convert count columns to numeric
numeric_cols1 <- c("T0", "D7", "PSU", "PSL")
rawcounts[numeric_cols1] <- lapply(rawcounts[numeric_cols1], function(x) as.numeric(as.character(x)))
#how many guides exist for each gene in the library
guides_per_gene <- rawcounts %>%
  group_by(Gene) %>%
  summarise(total_guides = n(), .groups = "drop")
#should probably save this cos it's interesting
#go back to counts data and treat anything <10 as 0
rawcounts_strict <- rawcounts %>%
  mutate(across(c(T0, D7, PSU, PSL), ~ ifelse(. < 10, 0, .)))
#how many guides drop out across all conditions now?
totaldropout <- rawcounts_strict %>%
  filter(T0 == 0 & D7 == 0 & PSU == 0 & PSL == 0)
#449 dropouts
missing_neuron_guides <- rawcounts_strict %>%
  filter(D7 == 0 & PSU == 0 & PSL == 0)
#2838 guides dropping out completely in neurons
#how many missing guides per gene (have no guides in neuronal populations)
missing_guides_per_gene <- missing_neuron_guides %>%
  group_by(Gene) %>%
  summarise(missing_guides = n(), .groups = "drop")
#combine to see which genes are no longer represented in neurons
gene_representation <- guides_per_gene %>%
  left_join(missing_guides_per_gene, by = "Gene") %>%
  mutate(missing_guides = ifelse(is.na(missing_guides), 0, missing_guides),
         guides_remaining = total_guides - missing_guides)
#genes with zero guides remaining in neuronal populations
genes_not_in_neurons <- gene_representation %>%
  filter(guides_remaining ==0)
#32 genes with no guides at all in neuronal populations

#now look at genes that aren't adequately represented in neuronal populations
genes_not_represented_in_neurons <- gene_representation %>%
  filter(guides_remaining <= 2)
#224 genes
#of these, 50 of them show up in the survival genes list
#load in survival stuff
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/")
survivalgenes = read.csv("droppedacross142128.csv")
#114 genes
#filter for more definite neuronal survival genes
survivalgenes_filtered <- survivalgenes %>%
  filter(mean_effect <= -2.5)
#110 genes
nosurvival <- intersect(
  genes_not_represented_in_neurons$Gene,
  survivalgenes_filtered$Gene
)
write.csv(nosurvival, "survivals_dropped.csv", row.names = F)
#of these, 50 of them show up in the survival genes list
gene_sets <- list(
  dropout = genes_not_represented_in_neurons$Gene,
  survivalgenes = survivalgenes_filtered$Gene)
venn.plot <- venn.diagram(
  x = gene_sets, filename = NULL,
  fill = c("#4DAF4A", "#984EA3"),
  alpha = 0.6,
  cex = 1.2,
  cat.cex = 1.2,
  cat.col = c("#4DAF4A", "#984EA3"))
grid.newpage()
grid.draw(venn.plot)



survival_gene_representation <- rawcounts_strict %>%
  filter(Gene %in% survivalgenes_filtered$Gene) %>%
  group_by(Gene) %>%
  summarise(
    total_guides = n(),
    guides_in_T0  = sum(T0  > 0),
    guides_in_D7  = sum(D7  > 0),
    guides_in_PSU = sum(PSU > 0),
    guides_in_PSL = sum(PSL > 0),
    .groups = "drop"
  ) %>%
  left_join(
    survivalgenes_filtered %>% select(Gene, mean_effect),
    by = "Gene"
  ) %>%
  arrange(guides_in_D7, guides_in_PSU, guides_in_PSL, mean_effect)
unsure = setdiff(survival_gene_representation$Gene, genes_not_represented_in_neurons$Gene)
#59 genes (unsure how we lost one)
write.csv(unsure, "survivals_notdropped.csv", row.names = F)

