rm (list = ls())
library(dplyr)
library(VennDiagram)
library(tidyverse)
library(ggplot2)
library(stringr)
library(ggrepel)
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/")
survivalgenes <- read.csv("droppedacross142128.csv")
setwd("/home/stoca738/Hugheslab_project/Caroline_Mageck/H1/mageck_nest/magecknest_test/")
dropouts <- read.csv("genestoremove.csv")
length(intersect(dropouts$Gene, survivalgenes$Gene))
#40
core_ess_gene_symbol=c('AAMP','AARS','AASDHPPT','ABCB7','ABCF1','ABT1','ACO2','ACTL6A','ACTR10','ACTR2','AHCY','ALDOA','ALG1','ALG11','ALG2','ANAPC10','ANAPC11','ANAPC4','ANAPC5','AP2S1','ARCN1','ARL2','ARMC7','ARPC4','ATP2A2','ATP5A1','ATP5B','ATP5F1','ATP5J2','ATP6AP1','ATP6V0B','ATP6V0C','ATP6V0D1','ATP6V1A','AURKA','BANF1','BCAS2','BIRC5','BMS1','BRD4','BRF1','BUB3','BUD31','BYSL','C21orf59','C3orf17','C9orf114','CAPZB','CARS','CASC5','CCDC115','CCDC84','CCNA2','CCT2','CCT3','CCT4','CCT5','CCT6A','CCT8','CDC123','CDC16','CDC20','CDC23','CDC27','CDC37','CDC45','CDC5L','CDC7','CDC73','CDIPT','CDK1','CDK9','CDT1','CEBPZ','CENPM','CENPN','CHAF1B','CHEK1','CHERP','CHMP2A','CHMP6','CIAO1','CINP','CIRH1A','CKAP5','CMPK1','CNOT3','COASY','COPA','COPB1','COPB2','COPE','COPS3','COPS6','CPSF2','CPSF3','CRCP','CRNKL1','CSE1L','CSNK1A1','CSNK2B','CSTF1','CSTF3','CTCF','CTDP1','CTU2',
                       'DAD1','DAP3','DARS','DBR1','DCTN5','DDB1','DDOST','DDX10','DDX18','DDX24','DDX27','DDX46','DDX47','DDX49','DDX52','DDX55','DDX56','DDX59','DDX6','DHDDS','DHX15','DHX36','DHX8','DHX9','DIEXF','DIS3','DKC1','DMAP1','DNAJA3','DNM2','DONSON','DTL','DTYMK','DYNC1H1','DYNLRB1','EBNA1BP2','ECD','ECT2','EEF2','EIF1','EIF2B1','EIF2B2','EIF2B3','EIF2B4','EIF2B5','EIF2S1','EIF2S2','EIF2S3','EIF3A','EIF3D','EIF3I','EIF4A1','EIF4A3','EIF5','EIF5B','EIF6','ELL','ELP5','ELP6','EPRS','ERCC3','ERH','ETF1','EXOSC2','EXOSC3','EXOSC5','EXOSC8','EXOSC9','FAM96B','FARSA','FARSB','FBL','FDX1L','FEN1','FNTB','FTSJ3','GAPDH','GARS','GEMIN5','GFER','GFM1','GGPS1','GINS1','GINS2','GNB1L','GNB2L1','GNL2','GNL3','GPKOW','GPN2','GPN3','GRPEL1','GRWD1','GSPT1','GTF2A2','GTF2B','GTF2H3','GTF3C1','GTF3C3','GTF3C4','GTF3C5','GTPBP4','GUK1','HARS','HAUS1','HAUS3','HAUS4','HAUS5','HAUS7','HCFC1','HEATR1','HINFP','HJURP','HMGCS1','HNRNPC','HNRNPK','HNRNPL','HNRNPU','HSD17B10','HSPA5','HSPA8','HSPA9','IARS','ICE1','IGBP1','IKBKAP','IMP4','INCENP','INTS2','INTS3','ISG20L2','KANSL3','KARS','KAT8','KIF11','KPNB1','KRR1','LARS','LAS1L','LRR1','LSM11','LSM2','LSM3','LSM4','LSM8','LTV1','MARS','MASTL','MAT2A','MAX','MCM2','MCM3','MCM5','MCM6','MCM7','MDN1','MED14','MED17','MED20','MED27','MED6','MED7','MED8','MED9','METTL14','METTL16','MFAP1','MIS12','MOB4','MRPL34','MRPL39','MRPL45','MRPS10','MRPS12','MRPS14','MRPS6','MTG2','MTOR','MYBBP1A','MZT1','NAA10','NAA15','NAA20','NAA25','NAA50','NAPA','NARFL','NARS','NCAPG','NCAPH','NCBP1','NCBP2','NCL','NDC80','NDUFAB1','NEDD1','NEDD8','NFS1','NGDN','NHP2L1','NIFK','NIP7','NKAP','NLE1','NMD3','NOC4L','NOL10','NOL12','NOL6','NOL8','NOL9','NOM1','NOP10','NOP14','NOP16','NOP2','NOP56','NOP58','NPAT','NPLOC4','NRF1','NSA2','NSMCE1','NUDT21','NUF2','NUP133','NUP155','NUP160','NUP85','NUP93','NUTF2',
                       'ORAOV1','ORC1','ORC6','OSGEP','PABPC1','PAK1IP1','PCF11','PCID2','PCNA','PDCD11','PDRG1','PES1','PFDN2','PFDN6','PHB','PLK1','PMPCA','PMPCB','PNO1','POLD1','POLD2','POLD3','POLE','POLE2','POLR1B','POLR1C','POLR2A','POLR2B','POLR2C','POLR2D','POLR2E','POLR2F','POLR2I','POLR2L','POLR3A','POLR3B','POLR3C','POLR3E','POLR3H','POLR3K','POP7','PPIL2','PPIL4','PPP1CB','PPP1R10','PPP4C','PPWD1','PRC1','PRIM1','PRMT1','PRPF19','PRPF3','PRPF31','PRPF38A','PRPF38B','PRPF4','PRPF6','PRPF8','PSMA1','PSMA2','PSMA3','PSMA5','PSMA7','PSMB1','PSMB2','PSMB3','PSMB4','PSMB7','PSMC1','PSMC3','PSMC4','PSMC5','PSMC6','PSMD1','PSMD12','PSMD14','PSMD2','PSMD3','PSMD4','PSMD6','PSMD7','PSMD8','PSMG4','PTPMT1','PWP2','QARS','RABGGTB','RAD51','RAE1','RAN','RANGAP1','RARS2','RBBP4','RBBP5','RBBP6','RBM14','RBM25','RBM8A','RBMX','RBX1','RCC1','RFC2','RFC3','RFC4','RFC5','RIOK1','RNGTT','RNPS1','RPA1','RPA2','RPF1','RPF2','RPL11','RPL12','RPL13','RPL13A','RPL14','RPL15','RPL17','RPL18','RPL18A','RPL19','RPL23','RPL26','RPL27','RPL27A','RPL29','RPL3','RPL31','RPL32','RPL34','RPL35','RPL37','RPL37A','RPL38','RPL4','RPL5','RPL6','RPL7A','RPL7L1','RPL9','RPLP0','RPLP2','RPN1','RPP14','RPP38','RPS10','RPS11','RPS12','RPS13','RPS14','RPS15','RPS15A','RPS16','RPS18','RPS19','RPS21','RPS23','RPS26','RPS27A','RPS29','RPS3','RPS4X','RPS8','RPS9','RPTOR','RRM1','RRM2','RRP1','RRP9','RSL1D1','RSL24D1','RTCB','RTFDC1','RUVBL1','RUVBL2','SACM1L','SAP18','SARS','SARS2','SART3','SDAD1','SEC61A1','SEH1L','SF1','SF3A2','SF3B3','SF3B6','SFPQ','SHFM1','SHQ1','SKIV2L2','SKP1','SLC35B1','SLC7A6OS','SLMO2','SMC1A','SMC2','SMC3','SMC4','SMC5','SMC6','SMNDC1','SMU1','SNAPC3','SNRNP200','SNRNP35','SNRNP70','SNRPA1','SNRPB','SNRPC','SNRPD1','SNRPD2','SNRPD3','SNRPF','SNW1','SOD1','SON','SPATA5','SPATA5L1','SPC24','SPC25','SPRTN','SRP72','SRRT','SRSF1','SRSF2','SRSF3','SRSF7','SSRP1','SSU72','SUMO2','SUPT16H','SUPV3L1','SYMPK','SYS1',
                       'TAF1B','TAF1C','TAMM41','TARDBP','TARS2','TBCD','TCP1','THOC5','TICRR','TIMM10','TIMM23','TINF2','TMEM199','TNPO3','TOMM20','TOMM22','TOMM40','TOP2A','TP53RK','TPI1','TPT1','TPX2','TRAPPC1','TRAPPC3','TRAPPC4','TRAPPC8','TRMT5','TSEN54','TSG101','TSR1','TSR2','TTC27','TTF1','TUBA1B','TUBB','TUBG1','TUBGCP3','TUBGCP4','TUT1','TWISTNB','TXNL4A','U2AF1','U2AF2','UBA1','UBA2','UBA3','UBA52','UBE2I','UBE2L3','UBL5','UFD1L','UPF1','UQCRC1','URB1','USP39','USP5','USPL1','UTP15','UTP20','UXT','VARS','VCP','VMP1','VPS13D','VPS25','WDHD1','WDR12','WDR18','WDR25','WDR3','WDR33','WDR43','WDR46','WDR55','WDR70','WDR74','WDR75','WDR77','WDR82','WDR92','WEE1','XAB2','XPO5','XRCC6','YARS','YARS2','YEATS4','YKT6','ZMAT5','ZNF131','ZPR1')
length(intersect(dropouts$Gene, core_ess_gene_symbol))
#18
length(intersect(survivalgenes$Gene, core_ess_gene_symbol))
#24
outliers <- read.delim("../noLL/noLL_suggested_outliers_by_gene.txt", stringsAsFactors = F)
nestnoLL <- read.delim("../noLL/noLL_PPI_False_outliers_removal_False.gene_summary.txt", stringsAsFactors = F)
nestsurvivalgenes <- nestnoLL %>%
  filter(Gene %in% survivalgenes$Gene)
nestsurvivalgenes <- nestsurvivalgenes %>%
  filter(!Gene %in% dropouts$Gene)
validnestsurvival <- nestsurvivalgenes %>%
  filter(!is.na(Pos.z) & is.finite(Pos.z),
         !is.na(CTRL.z) & is.finite(CTRL.z),
         sgRNA >= 3)

plot_df <- validnoLL %>%
  mutate(
    category = case_when(
      Gene %in% validnestsurvival$Gene ~ "Neuronal survival genes",
      TRUE ~ "Other"
    ))
ggplot(plot_df, aes(x = CTRL.beta, y = Pos.beta)) +
  geom_point(aes(color = category), alpha = 0.6) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "CTRL beta",
    y = "Pos beta (UL)",
    color = NULL,
    title = "Comparison of CTRL vs UL beta scores",
    subtitle = "Neuronal survival genes highlighted"
  ) +
  theme_classic()
ggplot(plot_df, aes(x = CTRL.beta, y = Pos.beta)) +
  geom_point(color = "grey80", alpha = 0.5) +
  geom_point(
    data = subset(plot_df, category == "Neuronal survival genes"),
    color = "red",
    size = 2
  ) +
  geom_text_repel(
    data = plot_df %>%
      filter(category == "Neuronal survival genes", Pos.beta > 0 | CTRL.beta > 0),
    aes(label = Gene),
    size = 3,
    max.overlaps = 50
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    x = "CTRL beta",
    y = "Pos beta (UL)",
    title = "Neuronal survival genes show distinct depletion patterns"
  ) +
  theme_classic()
#how many pos hits have outliers?
positive_hits <- plot_df %>%
  filter(Pos.beta > 0 | CTRL.beta > 0)
positive_hits_outliers <- positive_hits %>%
  left_join(outliers, by = "Gene")
positive_survival_outliers <- positive_hits_outliers %>%
  filter(category == "Neuronal survival genes")
outlier_summary <- outliers %>%
  group_by(Gene) %>%
  summarise(
    n_outlier_guides = n(),
    has_outlier = n_outlier_guides > 0,
    .groups = "drop"
  )
plot_df2 <- plot_df %>%
  left_join(outlier_summary, by = "Gene") %>%
  mutate(
    n_outlier_guides = ifelse(is.na(n_outlier_guides), 0, n_outlier_guides),
    has_outlier = ifelse(is.na(has_outlier), FALSE, has_outlier)
  )
ggplot(plot_df2, aes(x = CTRL.beta, y = Pos.beta)) +
  geom_point(color = "grey80", alpha = 0.5) +
  geom_point(
    data = plot_df2 %>%
      filter(category == "Neuronal survival genes", has_outlier),
    color = "blue",
    size = 2
  ) +
  geom_point(
    data = plot_df2 %>%
      filter(category == "Neuronal survival genes", !has_outlier),
    color = "red",
    size = 2
  ) +
  geom_text_repel(
    data = plot_df2 %>%
      filter(
        category == "Neuronal survival genes",
        !has_outlier,
        Pos.beta > 0 | CTRL.beta > 0
      ),
    aes(label = Gene),
    size = 3,
    max.overlaps = 50
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(
    title = "Neuronal survival genes show distinct depletion patterns",
    subtitle = "Blue = genes with sgRNA outliers",
    x = "D7 beta",
    y = "UL beta"
  ) +
  theme_classic()

validnoLL <- nestnoLL %>%
  filter(!is.na(Pos.z) & is.finite(Pos.z),
         !is.na(CTRL.z) & is.finite(CTRL.z),
         sgRNA >= 3)

ess_df <- validnoLL %>%
  mutate(
    is_essential = ifelse(Gene %in% survivalgenes$Gene, "Essential", "Other"),
    depletion = CTRL.beta)

ggplot(ess_df, aes(x = CTRL.beta, fill = is_essential)) +
  geom_density(alpha = 0.4) +
  labs(
    x = "Beta score",
    y = "Density",
    fill = NULL,
    title = "Essential genes show stronger depletion"
  ) +
  theme_classic()

ess_d7 <- validnoLL %>%
  mutate(
    is_essential = ifelse(Gene %in% survivalgenes$Gene, "Essential", "Other"),
    depletion = CTRL.beta)

ess_UL <- validnoLL %>%
  mutate(
    is_essential = ifelse(Gene %in% survivalgenes$Gene, "Essential", "Other"),
    depletion = Pos.beta)

ess_combined <- bind_rows(ess_UL, ess_d7)
ggplot(ess_combined, aes(x = is_essential, y = depletion, fill = is_essential)) +
  geom_boxplot() +
  labs(
    x = NULL,
    y = "Beta score",
    title = "Neuronal survival genes are depleted"
  ) +
  theme_classic()
