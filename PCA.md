# PCA of all samples (source material, in vitro skeletal muscle samples and reference skeletal muscle samples)

```R
#set working directory
setwd("/staging/leuven/stg_00134/GENERAL/rawdata/STARoutput/Ranalysis/")

#Libraries
library(tidyr)
library(magrittr)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(edgeR)
library(readr)
library(limma)
library(tidyverse)
library(limma)
library(RColorBrewer)
library(heatmaply) 
library(d3heatmap)

# read in count matric
counts = "/staging/leuven/stg_00134/GENERAL/rawdata/STARoutput/Ranalysis/Comparisons/allFC.counts.NOV.tsv"
counts <- as.matrix(read.csv(counts, sep="\t", row.names="X"))
colnames(counts) <- sub("^X", "", colnames(counts)) #to remove the Xs in the column names (this happens bc the column names start with a number)
dim(counts)
head(counts)

counts <- as.data.frame(counts)


# select all samples of interest
MT1 <- dplyr::select(counts, contains(c("Gene_ID","1D_MP","9D_MP","10D_MP","12D_MP_L9","12D_MP_L10","12D_MP_L11","12D_MP_L12","13D_MP","17D_MP","15D_MP","6D_MP13","6D_MP14","6D_MP15","6D_MP16","6D_MP17","6D_MP18","12D_MP_L13","12D_MP_L14","12D_MB_L1","12D_MB_L2","12D_MB_S1","D_MT","T_MP","T_MB","T_MT","E_MP","10F_FSM","13F_MP","33A","13A","20A_MB","46A_MB","47A_MB","49A_MB","50A_MB","52A_MB","54A_MB","60A_MB","43C_MB","39C","41C","45C","53C","24A_MB","55C","58C","47A","40C_MT","hPSC","T_FB","34A_T1","34A_T2", "61A", "62F", "63_3D")))
MT1.m <- as.matrix(MT1)
MT1

# get counts per million for PCA
DGE_MT1 <- DGEList(MT1.m)
cpmMT1 <- cpm(DGE_MT1)
log2.cpm_MT1 <- cpm(DGE_MT1, log=TRUE)
table(rowSums(DGE_MT1$counts==0)==418)
 

# filter out lowly expressed genes
keepersMT1 <- rowSums(cpmMT1>0)>=1
DGE_MT1.filt <- DGE_MT1[keepersMT1,]
dim(DGE_MT1.filt)

# normalize counts
log2.cpm_MT1 <- cpm(DGE_MT1.filt, log=TRUE)
DGE_MT1.filt.norm <- calcNormFactors(DGE_MT1.filt, method="TMM")
log2.cpm_MT1.filt.norm <- cpm(DGE_MT1.filt.norm, log=TRUE)

# factor the samples per cell type
groupMT1 <- factor(c("DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMP","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMP","TMT","TMT","TMT","TMT","TMT","TMT","TMT","TMT","TMT","TMT","TMT","TMT","TMT","TMT","TMT","FE","FE","FE","FE","FE","FE","FE","FE","FE","FE","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","Biopsy","SAT","SAT","SAT","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","PMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","CMB","PMT","PMT","CMT","CMT","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","FB","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF","MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", "MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"MF", 	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F",	"F", "3D","3D","3D","3D","3D","3D","3D"))

# check length of grouping variable to be same as number of samples
length(groupMT1)

# Add a sample label vector for labeling
samplelabels1 <- colnames(log2.cpm_MT1.filt.norm)

#hierarchical clustering: distance matrix 
distance <- dist(t(log2.cpm_MT1.filt.norm), method = "euclidean") 
clusters <- hclust(distance, method = "average")
plot(clusters, labels=samplelabels1, cex = 0.3)
```
![Hierarchical cluster](/Images/output_16_0.png)

```R    
#Principal component analysis of cpm filtered
pca.res <- prcomp(t(log2.cpm_MT1.filt.norm), scale.=F, retx=T)
summary(pca.res)
screeplot(pca.res)

pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per
pca.res.df <- as_tibble(pca.res$x)
```
![PCA RES](/Images/output_17_2.png)

```R  
#PCA without labels
ggplot(pca.res.df) +
  aes(x=PC1, y=PC2, color=groupMT1) +
  geom_point(size=2) +
  xlab(paste0("PC1 (",pc.per[1],"%",")")) + 
  ylab(paste0("PC2 (",pc.per[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_classic()
```
![PCA 1](/Images/output_18_0.png)

```R  
#PCA with labels
ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, color = groupMT1, label = samplelabels1) +
  geom_point(size = 2) +
  geom_text(vjust = -0.5, size = 3) +  # Add labels
  xlab(paste0("PC1 (", pc.per[1], "%", ")")) +
  ylab(paste0("PC2 (", pc.per[2], "%", ")")) +
  labs(title = "PCA plot",
       caption = paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_classic()
```
![PCA 2](/Images/output_19_0.png)

