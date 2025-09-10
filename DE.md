# Median filtering and generation of DE tables
###### This notebook shows the pipeline we have applied to each compariosn. In this case adult myogenic progenitors and immortalized myoblast are compared. Per comparison, only the included samples in line 43 need to be changes (sample information can be found in the supplementary excel). And the grouping of samples per category on line 71. 

```R
setwd("/staging/leuven/stg_00134/GENERAL/rawdata/STARoutput/Ranalysis/")

#Libraries
library(tximport)
library(tidyr)
library(magrittr)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(edgeR)
library(ensembldb)
library(rhdf5)
library(beepr)
library(EnsDb.Hsapiens.v86)
library(readr)
library(plotly)
library(limma)
library(gt)
library(DT)
library("VennDiagram")
library(tidyverse)
library(limma)
library(RColorBrewer)
library(gplots)
library(heatmaply) 
library(d3heatmap)

# read in the master count table
counts = "/staging/leuven/stg_00134/GENERAL/rawdata/STARoutput/Ranalysis/Comparisons/Master_Count.tsv"
counts <- as.matrix(read.csv(counts, sep="\t", row.names="X")) 
colnames(counts) <- sub("^X", "", colnames(counts)) #to remove the Xs in the column names (this happens bc the column names start with a number)
dim(counts)
head(counts)

counts <- as.data.frame(counts)

# select samples based on comparison (see sup table)
MT <- dplyr::select(counts, contains(c("Gene_ID","20A_MB1","20A_MB2","20A_MB3","20A_MB4","46A_MB1","46A_MB2","46A_MB3","47A_MB1","47A_MB2","49A_MB1","49A_MB2","49A_MB3","49A_MB4","49A_MB5","49A_MB6","49A_MB7","49A_MB8","50A_MB1","50A_MB2","52A_MB1","52A_MB2","52A_MB3","52A_MB4","52A_MB5","52A_MB6","52A_MB7","52A_MB8","52A_MB9","52A_MB10","54A_MB1","54A_MB2","60A_MB1","60A_MB2","43C_MB1","43C_MB2","13A_ASM1","13A_ASM2","13A_ASM3")))
MT

# save the column names as labels 
cols <- colnames(MT)
cols

# convert to matrix 
MT.m <- as.matrix(MT)
MT.m

# differential gene expression analysis start
DGE_MT <- DGEList(MT.m)
cpmMT <- cpm(DGE_MT)
log2.cpm_MT <- cpm(DGE_MT, log=TRUE)
table(rowSums(DGE_MT$counts==0)==38) # change the seocnd number depending on the total samples you have

# filter loly expressed genes
keepersMT <- rowSums(cpmMT>0)>=1
DGE_MT.filt <- DGE_MT[keepersMT,]
dim(DGE_MT.filt)

# Normalization
log2.cpm_MT <- cpm(DGE_MT.filt, log=TRUE)
DGE_MT.filt.norm <- calcNormFactors(DGE_MT.filt, method="TMM")
log2.cpm_MT.filt.norm <- cpm(DGE_MT.filt.norm, log=TRUE)

# grouping variables to define 
groupMT <- factor(c("P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","P","A","A","A"))
length(groupMT)
groupMT

# save sample labels of the normalized and filtered matrix
samplelabels <- colnames(log2.cpm_MT.filt.norm)
samplelabels

# design matrix per group
MTdesign <- model.matrix(~0 + groupMT)
colnames(MTdesign) <- levels(groupMT)
v.DGE_MT.filt.norm <- voom(DGE_MT.filt.norm, MTdesign, plot = TRUE)
MTfit <- lmFit(v.DGE_MT.filt.norm, MTdesign)

# Decide on the direction of up- and downregulated genes
MTcontrast.matrix <- makeContrasts(Difference = P - A,
                                 levels=MTdesign)
MTfits <- contrasts.fit(MTfit, MTcontrast.matrix)
MTebFit <- eBayes(MTfits)

# number of significantly up and down regulated genes where p<0.05>
MTresults <- decideTests(MTebFit, method="global", adjust.method="BH", p.value=0.05)
summary(MTresults)
MTTopHits <- topTable(MTebFit, adjust ="BH", coef=1, number=12203, sort.by="P") # change number to total significant DEGs found on line 93

# TopHits table of logFC
MTTopHits.df <- MTTopHits %>% as_tibble(rownames = "Gene_ID")

# dataframe of DEGs
MT_DG <- v.DGE_MT.filt.norm$E[MTresults[,1] !=0,]
MT_DG.df <- as_tibble(MT_DG, rownames = "Gene_ID")

# Downregulation
# Define thresholds
fold_change_threshold <- 0  # Use a negative threshold for downregulation
p_value_threshold <- 0.05

# Subset data frame for downregulated genes
downregulatedTopHits <- MTTopHits.df[MTTopHits.df$logFC < fold_change_threshold & MTTopHits.df$adj.P.Val < p_value_threshold, ]
downregulatedTopHits %>% arrange((logFC))

# Downregulation
# Define thresholds
fold_change_threshold <- 0  # Use a negative threshold for downregulation
p_value_threshold <- 0.05

# Subset data frame for downregulated genes
upregulatedTopHits <- MTTopHits.df[MTTopHits.df$logFC > fold_change_threshold & MTTopHits.df$adj.P.Val < p_value_threshold, ]
upregulatedTopHits %>% arrange(desc(logFC))
```


# median filtering


```R
# change this dependent on the comparison
TMT_cols <- c("20A_MB1","20A_MB2","20A_MB3","20A_MB4","46A_MB1","46A_MB2","46A_MB3","47A_MB1","47A_MB2","49A_MB1","49A_MB2","49A_MB3","49A_MB4","49A_MB5","49A_MB6","49A_MB7","49A_MB8","50A_MB1","50A_MB2","52A_MB1","52A_MB2","52A_MB3","52A_MB4","52A_MB5","52A_MB6","52A_MB7","52A_MB8","52A_MB9","52A_MB10","54A_MB1","54A_MB2","60A_MB1","60A_MB2","43C_MB1","43C_MB2")
MF_cols <- c("13A_ASM1","13A_ASM2","13A_ASM3")
length(TMT_cols)
length(MF_cols)

cpmMT.df <- as.data.frame(cpmMT)
cpmMT.df <- rownames_to_column(cpmMT.df, var="Gene_ID")

cpmMT.df #ok

cpmMTmed.df <- cpmMT.df %>%
  rowwise() %>%
  mutate(
    median_TMT = median(c_across(all_of(TMT_cols))),
    median_MF = median(c_across(all_of(MF_cols)))

cpmMTmed.df #ok

# all tried filters
medTMTmore1 <- cpmMTmed.df %>%
  filter(median_TMT > 1)
medTMTmore1

medTMTless1 <- cpmMTmed.df %>%
  filter(median_TMT < 1)
medTMTless1

medMFless1 <- cpmMTmed.df %>%
  filter(median_MF < 1)
medMFless1

medMFmore1 <- cpmMTmed.df %>%
  filter(median_MF > 1)
medMFmore1

medMFmore5 <- cpmMTmed.df %>%
  filter(median_MF > 5)
medMFmore5

medTMTmore5 <- cpmMTmed.df %>%
  filter(median_TMT > 5)
medTMTmore5

medMFless5 <- cpmMTmed.df %>%
  filter(median_MF < 5)
medMFless5

medTMTless5 <- cpmMTmed.df %>%
  filter(median_TMT < 5)
medTMTless5
```


# filter 


```R
# downregulatedTopHits
# upregulatedTopHits

# downregulated genes so - in TMT less1 only one filter 
DRtophitsFILTER1 <- inner_join(downregulatedTopHits, medTMTless1, by = "Gene_ID")
DRtophitsFILTER1

DRtophitsFILTER2 <- inner_join(DRtophitsFILTER1,medMFmore1, by = "Gene_ID")
DRtophitsFILTER2

DRtophitsFILTER3 <- inner_join(DRtophitsFILTER1,medMFmore5, by = "Gene_ID")
DRtophitsFILTER3

# downregulated genes so - in TMT less1 only one filter 
#DRtophitsFILTER4 <- inner_join(downregulatedTopHits, medTMTless5, by = "Gene_ID")
#DRtophitsFILTER4

# downregulated genes so - in TMT less1 only one filter 
#DRtophitsFILTER5 <- inner_join(DRtophitsFILTER4, medMFmore5, by = "Gene_ID")
#DRtophitsFILTER5

# upregulated genes so - in MF less1 only one filter 
URtophitsFILTER1 <- inner_join(upregulatedTopHits, medMFless1, by = "Gene_ID")
URtophitsFILTER1 %>% arrange(desc(logFC))

# upregulated genes so - in MF less1 only one filter 
URtophitsFILTER2 <- inner_join(URtophitsFILTER1, medTMTmore1, by = "Gene_ID")
URtophitsFILTER2 %>% arrange(desc(logFC))

# upregulated genes so - in MF less1 only one filter 
URtophitsFILTER3 <- inner_join(URtophitsFILTER1, medTMTmore5, by = "Gene_ID")
URtophitsFILTER3 %>% arrange(desc(logFC))

# upregulated genes so - in MF less1 only one filter 
#URtophitsFILTER4 <- inner_join(upregulatedTopHits, medMFless5, by = "Gene_ID")
#URtophitsFILTER4 %>% arrange(desc(logFC))
```


```R
# upregulated genes so - in MF less1 only one filter 
#URtophitsFILTER5 <- inner_join(URtophitsFILTER4, medTMTmore5, by = "Gene_ID")
#URtophitsFILTER5 %>% arrange(desc(logFC))
```

# what I will work with


```R
DRtophitsFILTER3

URtophitsFILTER3


UPgenefilter <- URtophitsFILTER3$Gene_ID
DOWNgenefilter <- DRtophitsFILTER3$Gene_ID


upTH <- upregulatedTopHits[upregulatedTopHits$Gene_ID %in% UPgenefilter, ]
upTH

downTH <- downregulatedTopHits[downregulatedTopHits$Gene_ID %in% DOWNgenefilter, ]
downTH

fullTH <- bind_rows(upTH,downTH)
fullTH %>% arrange(desc(logFC))

fullTHcpm <- inner_join(fullTH, cpmMT.df, by = "Gene_ID")
upTHcpm <- inner_join(upTH, cpmMT.df, by = "Gene_ID")
downTHcpm <- inner_join(downTH, cpmMT.df, by = "Gene_ID")
fullTHcpm

write.table(fullTHcpm, "/staging/leuven/stg_00134/Margaux/0710/IC5/TDMTAfullTHcpm.tsv", sep="\t", col.names=NA, quote=FALSE)
write.table(upTHcpm, "/staging/leuven/stg_00134/Margaux/0710/IC5/TDMTAupTHcpm.tsv", sep="\t", col.names=NA, quote=FALSE)
write.table(downTHcpm, "/staging/leuven/stg_00134/Margaux/0710/IC5/TDMTAdownTHcpm.tsv", sep="\t", col.names=NA, quote=FALSE)
```

# volcanoplot


```R
DGEvplot2 <- ggplot(MTTopHits.df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", Gene_ID)) +
  geom_point(aes(color = Gene_ID %in% c("DKK1", "PTX3", "RGS4", "GALNT5", "ADAMTSL1", "CA12", "RAB3B", "ALPL", "SEMA7A", "FGF5", "LYPD6B", "NT5E", "ASPM", "DEPDC1","DLK1", "FOSB", "SPARCL1", "CXCL14", "EDN3", "LPL", "HSPA6", "CHRDL2", "HEY2", "CALCR", "ITIH5", "RGMA", "MYBPC1", "TEX101")), size=2) +
  scale_color_manual(values = c("TRUE" = "pink", "FALSE" = "black")) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", size=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", size=1) +
  geom_text(data = subset(MTTopHits.df, Gene_ID %in% c("DKK1", "PTX3", "RGS4", "GALNT5", "ADAMTSL1", "CA12", "RAB3B", "ALPL", "SEMA7A", "FGF5", "LYPD6B", "NT5E", "ASPM", "DEPDC1","DLK1", "FOSB", "SPARCL1", "CXCL14", "EDN3", "LPL", "HSPA6", "CHRDL2", "HEY2", "CALCR", "ITIH5", "RGMA", "MYBPC1", "TEX101")),
            aes(label = Gene_ID), vjust = -1, hjust = 1) +
  labs(title="Volcano plot DGE IMM-ABiop") +
  theme(legend.position = "none")

DGEvplot3 <- DGEvplot2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"))
```

# Gene set filtering 


```R
MYO_THup <- upTH %>% filter(Gene_ID %in% c("ABLIM1","ACHE","ACSL1","ACTA1","ACTC1","ACTN2","ACTN3","ADAM12","ADCY9","AEBP1","AGL","AGRN","AK1","AKT2","ANKRD2","APLNR","APOD","APP","ATP2A1","ATP6AP1","BAG1","BDKRB2","BHLHE40","BIN1","CACNA1H","CACNG1","CAMK2B","CASQ1","CASQ2","CAV3","CD36","CDH13","CDKN1A","CFD","CHRNA1","CHRNB1","CHRNG","CKB","CKM","CKMT2","CLU","CNN3","COL15A1","COL1A1","COL3A1","COL4A2","COL6A2","COL6A3","COX6A2","COX7A1","CRAT","CRYAB","CSRP3","CTF1","DAPK2","DES","DMD","DMPK","DTNA","EFS","EIF4A2","ENO3","EPHB3","ERBB3","FABP3","FDPS","FGF2","FHL1","FKBP1B","FLII","FOXO4","FST","FXYD1","GAA","GABARAPL2","GADD45B","GJA5","GNAO1","GPX3","GSN","HBEGF","HDAC5","HRC","HSPB2","HSPB8","IFRD1","IGF1","IGFBP3","IGFBP7","ITGA7","ITGB1","ITGB4","ITGB5","KCNH1","KCNH2","KIFC3","KLF5","LAMA2","LARGE1","LDB3","LPIN1","LSP1","MAPK12","MAPRE3","MB","MEF2A","MEF2C","MEF2D","MRAS","MYBPC3","MYBPH","MYF6","MYH1","MYH11","MYH2","MYH3","MYH4","MYH7","MYH8","MYH9","MYL1","MYL2","MYL3","MYL4","MYL6B","MYL7","MYLK","MYL11","MYO1C","MYOG","MYOM1","MYOM2","MYOZ1","NAV2","NCAM1","NOS1","NOTCH1","NQO1","OCEL1","PC","PDE4DIP","PDLIM7","PFKM","PICK1","PKIA","PLXNB2","PPFIA4","PPP1R3C","PRNP","PSEN2","PTGIS","PTP4A3","PVALB","PYGM","RB1","REEP1","RIT1","RYR1","SCD","SCHIP1","SGCA","SGCD","SGCG","SH2B1","SH3BGR","SIRT2","SLC6A8","SLN","SMTN","SOD3","SORBS1","SORBS3","SPARC","SPDEF","SPEG","SPHK1","SPTAN1","SSPN","DENND2B","STC2","SVIL","SYNGR2","TAGLN","TCAP","TEAD4","TGFB1","TNNC1","TNNC2","TNNI1","TNNI2","TNNT1","TNNT2","TNNT3","TPD52L1","TPM2","TPM3","TSC2","VIPR1","WWTR1","NEO1", "TCF3", "MAPK14", "CTNNA2", "MEF2B", "CDH15", "MYF5", "ABL1", "BNIP2", "CDH2", "CDH4", "TCF12", "CTNNA1", "SPAG9", "NTN3", "BOC", "TCF4", "CDON", "CTNNB1", "CDC42", "MYOD1", "MAP2K6", "MAPK11"))
MYO_THup


MYO_THdown <- downTH %>% filter(Gene_ID %in% c("ABLIM1","ACHE","ACSL1","ACTA1","ACTC1","ACTN2","ACTN3","ADAM12","ADCY9","AEBP1","AGL","AGRN","AK1","AKT2","ANKRD2","APLNR","APOD","APP","ATP2A1","ATP6AP1","BAG1","BDKRB2","BHLHE40","BIN1","CACNA1H","CACNG1","CAMK2B","CASQ1","CASQ2","CAV3","CD36","CDH13","CDKN1A","CFD","CHRNA1","CHRNB1","CHRNG","CKB","CKM","CKMT2","CLU","CNN3","COL15A1","COL1A1","COL3A1","COL4A2","COL6A2","COL6A3","COX6A2","COX7A1","CRAT","CRYAB","CSRP3","CTF1","DAPK2","DES","DMD","DMPK","DTNA","EFS","EIF4A2","ENO3","EPHB3","ERBB3","FABP3","FDPS","FGF2","FHL1","FKBP1B","FLII","FOXO4","FST","FXYD1","GAA","GABARAPL2","GADD45B","GJA5","GNAO1","GPX3","GSN","HBEGF","HDAC5","HRC","HSPB2","HSPB8","IFRD1","IGF1","IGFBP3","IGFBP7","ITGA7","ITGB1","ITGB4","ITGB5","KCNH1","KCNH2","KIFC3","KLF5","LAMA2","LARGE1","LDB3","LPIN1","LSP1","MAPK12","MAPRE3","MB","MEF2A","MEF2C","MEF2D","MRAS","MYBPC3","MYBPH","MYF6","MYH1","MYH11","MYH2","MYH3","MYH4","MYH7","MYH8","MYH9","MYL1","MYL2","MYL3","MYL4","MYL6B","MYL7","MYLK","MYL11","MYO1C","MYOG","MYOM1","MYOM2","MYOZ1","NAV2","NCAM1","NOS1","NOTCH1","NQO1","OCEL1","PC","PDE4DIP","PDLIM7","PFKM","PICK1","PKIA","PLXNB2","PPFIA4","PPP1R3C","PRNP","PSEN2","PTGIS","PTP4A3","PVALB","PYGM","RB1","REEP1","RIT1","RYR1","SCD","SCHIP1","SGCA","SGCD","SGCG","SH2B1","SH3BGR","SIRT2","SLC6A8","SLN","SMTN","SOD3","SORBS1","SORBS3","SPARC","SPDEF","SPEG","SPHK1","SPTAN1","SSPN","DENND2B","STC2","SVIL","SYNGR2","TAGLN","TCAP","TEAD4","TGFB1","TNNC1","TNNC2","TNNI1","TNNI2","TNNT1","TNNT2","TNNT3","TPD52L1","TPM2","TPM3","TSC2","VIPR1","WWTR1","NEO1", "TCF3", "MAPK14", "CTNNA2", "MEF2B", "CDH15", "MYF5", "ABL1", "BNIP2", "CDH2", "CDH4", "TCF12", "CTNNA1", "SPAG9", "NTN3", "BOC", "TCF4", "CDON", "CTNNB1", "CDC42", "MYOD1", "MAP2K6", "MAPK11"))
MYO_THdown

write.table(MYO_THdown, "/staging/leuven/stg_00134/Margaux/0710/IC5/MYO_THdown.tsv", sep="\t", col.names=NA, quote=FALSE)
write.table(MYO_THup, "/staging/leuven/stg_00134/Margaux/0710/IC5/MYO_THup.tsv", sep="\t", col.names=NA, quote=FALSE)

# human TF
humanTF <- read_csv("/staging/leuven/stg_00134/Margaux/HumanTFs.csv")
humanTF <- as.data.frame(humanTF)
colnames(humanTF) <- c("Gene_ID")
humanTF

TF_THdown <- inner_join(downTH, humanTF, by = "Gene_ID")
TF_THup <- inner_join(upTH, humanTF, by = "Gene_ID")

TF_THdown

TF_THup %>% arrange(desc(logFC))

# Epifactors
Epi <- read_csv("/staging/leuven/stg_00134/Margaux/EpiGenes_main.csv")
Epi <- Epi[, "HGNC_symbol"]
colnames(Epi) <- c("Gene_ID")
Epi

EPI_THdown <- inner_join(downTH, Epi, by = "Gene_ID")
EPI_THup <- inner_join(upTH, Epi, by = "Gene_ID")

EPI_THdown

EPI_THup

# Metabolism
Metabolism <- read_csv("/staging/leuven/stg_00134/Margaux/metabolism_gene_list.csv")
colnames(Metabolism) <- c("Gene_ID")
Metabolism

MET_THdown <- inner_join(downTH, Metabolism, by = "Gene_ID")
MET_THup <- inner_join(upTH, Metabolism, by = "Gene_ID")

MET_THdown

MET_THup

write.table(TF_THdown, "/staging/leuven/stg_00134/Margaux/0710/IC5/TF_THdown.tsv", sep="\t", col.names=NA, quote=FALSE)
write.table(TF_THup, "/staging/leuven/stg_00134/Margaux/0710/IC5/TF_THup.tsv", sep="\t", col.names=NA, quote=FALSE)
write.table(EPI_THdown, "/staging/leuven/stg_00134/Margaux/0710/IC5/EPI_THdown.tsv", sep="\t", col.names=NA, quote=FALSE)
write.table(EPI_THup, "/staging/leuven/stg_00134/Margaux/0710/IC5/EPI_THup.tsv", sep="\t", col.names=NA, quote=FALSE)
write.table(MET_THdown, "/staging/leuven/stg_00134/Margaux/0710/IC5/MET_THdown.tsv", sep="\t", col.names=NA, quote=FALSE)
write.table(MET_THup, "/staging/leuven/stg_00134/Margaux/0710/IC5/MET_THup.tsv", sep="\t", col.names=NA, quote=FALSE)
```

