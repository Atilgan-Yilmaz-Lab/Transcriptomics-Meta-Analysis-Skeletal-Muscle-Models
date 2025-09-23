### notebook looking at PVCA and PCA of hPSC and DMT


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

counts = "/staging/leuven/stg_00134/GENERAL/rawdata/STARoutput/Ranalysis/Comparisons/allFC.counts.NOV.tsv"
counts <- as.matrix(read.csv(counts, sep="\t", row.names="X"))
colnames(counts) <- sub("^X", "", colnames(counts)) #to remove the Xs in the column names (this happens bc the column names start with a number)
counts <- as.data.frame(counts)
dim(counts)
head(counts)

MT1 <- dplyr::select(counts, contains(c("Gene_ID","hPSC", "D_MT", "12D_MB_L1", "12D_MB_L2", "12D_MB_S1")))
MT1

colnames(MT1)


MT1.m <- as.matrix(MT1)

DGE_MT1 <- DGEList(MT1.m)
cpmMT1 <- cpm(DGE_MT1)
log2.cpm_MT1 <- cpm(DGE_MT1, log=TRUE)
table(rowSums(DGE_MT1$counts==0)==35)

keepersMT1 <- rowSums(cpmMT1>0)>=1
DGE_MT1.filt <- DGE_MT1[keepersMT1,]
dim(DGE_MT1.filt)

log2.cpm_MT1 <- cpm(DGE_MT1.filt, log=TRUE)
DGE_MT1.filt.norm <- calcNormFactors(DGE_MT1.filt, method="TMM")
log2.cpm_MT1.filt.norm <- cpm(DGE_MT1.filt.norm, log=TRUE)

groupMT1 <- factor(c("hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","hPSC","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT","DMT"))
length(groupMT1)

samplelabels1 <- colnames(log2.cpm_MT1.filt.norm)

#hierarchical clustering: distance matrix 
distance <- dist(t(log2.cpm_MT1.filt.norm), method = "euclidean") 
clusters <- hclust(distance, method = "average")
plot(clusters, labels=samplelabels1, cex = 0.6)
```
![Hierarchical cluster](/Images/PVCA/output_12_0.png)   
    
```R
#Principal component analysis of cpm filtered
pca.res <- prcomp(t(log2.cpm_MT1.filt.norm), scale.=F, retx=T)
summary(pca.res)
screeplot(pca.res)

pc.var<-pca.res$sdev^2
pc.per<-round(pc.var/sum(pc.var)*100, 1)
pc.per
pca.res.df <- as_tibble(pca.res$x)


# Extract the prefix (e.g., "1D", "13A", etc.) into a new column
prefix <- str_extract(samplelabels1, "^[0-9]+[A-Z]")
prefix


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
![Hierarchical cluster](/Images/PVCA/output_15_0.png) 

```R
D_PCA <- ggplot(pca.res.df) +
  aes(x = PC1, y = PC2, color = prefix, label = groupMT1) +
  geom_point(size = 2) +
  geom_text(vjust = -0.5, size = 3) +  # Add labels
  xlab(paste0("PC1 (", pc.per[1], "%", ")")) +
  ylab(paste0("PC2 (", pc.per[2], "%", ")")) +
  labs(title = "PCA plot",
       caption = paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_classic()

D_PCA
```
    
![Hierarchical cluster](/Images/PVCA/output_17_0.png) 
    

```R
#creating a PCA small multiples chart + we are using the %>% pipe function from the magrittr package for the first time
pca.res.df <- pca.res$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = samplelabels1)

#pivoting the data frame for tidy data
pca.pivot <- pivot_longer(pca.res.df, # dataframe to be pivoted
                          cols = PC1:PC4, # column names to be stored as a SINGLE variable
                          names_to = "PC", # name of that new variable (column)
                          values_to = "loadings") # name of the new variable (column) storing the data


PCvar <- ggplot(pca.pivot) +
  aes(x=sample, y=loadings) + # you could iteratively 'paint' different covariates onto this plot using the 'fill' aes
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()
PCvar
```
    
![Hierarchical cluster](/Images/PVCA/output_19_0.png) 
    

```R
library(variancePartition)

log2.cpm_MT1.filt.norm

expr_mat <- log2.cpm_MT1.filt.norm
stopifnot(is.matrix(expr_mat))

colnames(expr_mat)

sample_info <- read.csv("/staging/leuven/stg_00134/Margaux/metadata.csv")
sample_info <- sample_info %>% select(Sample, Lab, Substrate, CellType, PSCtype, Seq, Protocol)
sample_info

rownames(sample_info) <- sample_info$SampleID
sample_info$SampleID <- NULL

# Ensure covariates are factors (categorical)
sample_info$Lab <- factor(sample_info$Lab)
sample_info$Substrate <- factor(sample_info$Substrate)
sample_info$CellType <- factor(sample_info$CellType)
sample_info$PSCtype <- factor(sample_info$PSCtype)
sample_info$Seq <- factor(sample_info$Seq)
sample_info$Protocol <- factor(sample_info$Protocol)

rownames(sample_info) <- sample_info$Sample
sample_info$Sample <- NULL
sample_info

form <- ~ (1|Lab) + (1|Substrate) + (1|CellType) + (1|PSCtype) + (1|Seq) + (1|Protocol)

vp <- fitExtractVarPartModel(expr_mat, form, sample_info)
vp <- sortCols(vp)
plotVarPart(vp)

PVCAplot <- plotVarPart(vp)

PVCAplot
```
![Hierarchical cluster](/Images/PVCA/output_33_0.png) 
