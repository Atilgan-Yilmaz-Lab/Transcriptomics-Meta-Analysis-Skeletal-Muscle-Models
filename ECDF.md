
### This workbook calculates the percentages of genes passing the appropriate cpm filters cpm >5 and cpm <1. Apply both filters separately to gneerate 2 combined ECDF plots which can be found in supplementary figure 1D and E
```R
# set working directory
setwd("/staging/leuven/stg_00134/Margaux/Revision/BatchCorrection/D1")

# load libraries
library(tidyverse)

# Read in comparison tables of interest
DMP_AMP <- read_csv("DMP_AMP.csv")
DMP_FMP <- read_csv("DMP_FMP.csv")
DMT_AMF <- read_csv("DMT_AMF.csv")
DMT_FETB <- read_csv("DMT_FETB.csv")
TDMB_AMP <- read_csv("TDMB_AMP.csv")
TDMB_FMP <- read_csv("TDMB_FMP.csv")
TDMT_AMF <- read_csv("TDMT_AMF.csv")
TDMT_FETB <- read_csv("TDMT_FETB.csv")

# set cpm threshold which was used for downregulated gene filtering 
threshold <- 1
# threshold <- 5 

# filter all comparison count tables for downregulated genes in the in vitro models
# Calculates the percentage of the samples for which a gene is passing the cpm filter
# swap to threshold <- 5 filter for the upregulated equivalent 
1DMP_AMP_pct <- DMP_AMP %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("D")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

1DMP_FMP_pct <- DMP_FMP %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("D")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

1DMT_AMF_pct <- DMT_AMF %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("D")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

1DMT_FETB_pct <- DMT_FETB %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("D")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

1TDMB_AMP_pct <- TDMB_AMP %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("T")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

1TDMB_FMP_pct <- TDMB_FMP %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("T")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

1TDMT_AMF_pct <- TDMT_AMF %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("T")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

1TDMT_FETB_pct <- TDMT_FETB %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("T")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

  # example of code for the upregukated cpm >5 filter
#1TDMT_FETB_pct <- TDMT_FETB %>%
#  rowwise() %>%
#  mutate(PctSamplesPassing = mean(c_across(contains("T")) >= threshold) * 100) %>%
#  ungroup() %>%
#  select(gene, PctSamplesPassing)  


# check dimensions of the 
dim(1DMP_AMP_pct)
dim(1DMP_FMP_pct)
dim(1DMT_AMF_pct)
dim(1DMT_FETB_pct)
dim(1TDMB_AMP_pct)
dim(1TDMB_FMP_pct)
dim(1TDMT_AMF_pct)
dim(1TDMT_FETB_pct)

# Percentage of genes where more than 80% of samples pass the filter
more80DMP_AMP <- 1DMP_AMP_pct %>% filter(PctSamplesPassing >80)
dim(more80DMP_AMP)

more80DMP_FMP <- 1DMP_FMP_pct %>% filter(PctSamplesPassing >80)
dim(more80DMP_FMP)

more80DMT_AMF <- 1DMT_AMF_pct %>% filter(PctSamplesPassing >80)
dim(more80DMT_AMF)

more80DMT_FETB <- 1DMT_FETB_pct %>% filter(PctSamplesPassing >80)
dim(more80DMT_FETB)

more80TDMB_AMP <- 1TDMB_AMP_pct %>% filter(PctSamplesPassing >80)
dim(more80TDMB_AMP)

more80TDMB_FMP <- 1TDMB_FMP_pct %>% filter(PctSamplesPassing >80)
dim(more80TDMB_FMP)

more80TDMT_AMF <- 1TDMT_AMF_pct %>% filter(PctSamplesPassing >80)
dim(more80TDMT_AMF)

more80TDMT_FETB <- 1TDMT_FETB_pct %>% filter(PctSamplesPassing >80)
dim(more80TDMT_FETB)


# Percentage of genes where more than 60% of samples pass the filter
more60DMP_AMP <- 1DMP_AMP_pct %>% filter(PctSamplesPassing >60)
dim(more60DMP_AMP)

more60DMP_FMP <- 1DMP_FMP_pct %>% filter(PctSamplesPassing >60)
dim(more60DMP_FMP)

more60DMT_AMF <- 1DMT_AMF_pct %>% filter(PctSamplesPassing >60)
dim(more60DMT_AMF)

more60DMT_FETB <- 1DMT_FETB_pct %>% filter(PctSamplesPassing >60)
dim(more60DMT_FETB)

more60TDMB_AMP <- 1TDMB_AMP_pct %>% filter(PctSamplesPassing >60)
dim(more60TDMB_AMP)

more60TDMB_FMP <- 1TDMB_FMP_pct %>% filter(PctSamplesPassing >60)
dim(more60TDMB_FMP)

more60TDMT_AMF <- 1TDMT_AMF_pct %>% filter(PctSamplesPassing >60)
dim(more60TDMT_AMF)

more60TDMT_FETB <- 1TDMT_FETB_pct %>% filter(PctSamplesPassing >60)
dim(more60TDMT_FETB)


# Median percentage of samples per gene passing the filter
medianDMP_AMP <- 1DMP_AMP_pct %>% summarise(median(PctSamplesPassing))
medianDMP_AMP

medianDMP_FMP <- 1DMP_FMP_pct %>% summarise(median(PctSamplesPassing))
medianDMP_FMP

medianDMT_AMF <- 1DMT_AMF_pct %>% summarise(median(PctSamplesPassing))
medianDMT_AMF

medianDMT_FETB <- 1DMT_FETB_pct %>% summarise(median(PctSamplesPassing))
medianDMT_FETB

medianTDMB_AMP <- 1TDMB_AMP_pct %>% summarise(median(PctSamplesPassing))
medianTDMB_AMP

medianTDMB_FMP <- 1TDMB_FMP_pct %>% summarise(median(PctSamplesPassing))
medianTDMB_FMP

medianTDMT_AMF <- 1TDMT_AMF_pct %>% summarise(median(PctSamplesPassing))
medianTDMT_AMF

medianTDMT_FETB <- 1TDMT_FETB_pct %>% summarise(median(PctSamplesPassing))
medianTDMT_FETB

# individual density plots
DMP_AMPdensity <- ggplot(1DMP_AMP_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (DMP_AMP)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

DMP_FMPdensity <- ggplot(1DMP_FMP_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (DMP_FMP)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

DMT_AMFdensity <- ggplot(1DMT_AMF_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (DMT_AMF)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

DMT_FETBdensity <- ggplot(1DMT_FETB_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (DMT_FETB)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

TDMB_AMPdensity <- ggplot(1TDMB_AMP_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (TDMB_AMP)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

TDMB_FMPdensity <- ggplot(1TDMB_FMP_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (TDMB_FMP)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

TDMT_AMFdensity <- ggplot(1TDMT_AMF_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (TDMT_AMF)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

TDMT_FETBdensity <- ggplot(1TDMT_FETB_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (TDMT_FETB)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()



# individual ECDF plots
DMP_AMPecdf <- ggplot(1DMP_AMP_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (DMP_AMP)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

DMP_FMPecdf <- ggplot(1DMP_FMP_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (DMP_FMP)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

DMT_AMFecdf <- ggplot(1DMT_AMF_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (DMT_AMF)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

DMT_FETBecdf <- ggplot(1DMT_FETB_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (DMT_FETB)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

TDMB_AMPecdf <- ggplot(1TDMB_AMP_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (TDMB_AMP)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

TDMB_FMPecdf <- ggplot(1TDMB_FMP_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (TDMB_FMP)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

TDMT_AMFecdf <- ggplot(1TDMT_AMF_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (TDMT_AMF)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

TDMT_FETBecdf <- ggplot(1TDMT_FETB_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (TDMT_FETB)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()



# Combine the individually filtered tables and investigate the total numbe rof genes passing the cpm filter in the in vitro direction
1Combined_pct <- bind_rows(DMP_AMP_pct, DMP_FMP_pct, DMT_AMF_pct, DMT_FETB_pct, TDMB_AMP_pct, TDMB_FMP_pct, TDMT_AMF_pct, TDMT_FETB_pct)
dim(1Combined_pct)
head(1Combined_pct)

# save this combined table as a csv
write.csv(1Combined_pct, "Combined_pctD1.csv")

# Number of genes passing for more than 80% of samples in the in vitro conditions
more80 <- 1Combined_pct %>% filter(PctSamplesPassing >80)
more80


# Number of genes passing for more than 60% of samples in the in vitro conditions
more60 <- 1Combined_pct %>% filter(PctSamplesPassing >60)
more60

# Median percenatge of samples per gene passing the filter
medianPCtpassing <- 1Combined_pct %>% summarise(median(PctSamplesPassing))
medianPCtpassing

# density plot for combined table
1Combined_density <- ggplot(1Combined_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "pink", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()
1Combined_density

# combined table ECDF
1Combined_ecdf <- ggplot(1Combined_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()
1Combined_ecdf

# Repeat the cpm <1 threshold but this time for the in vivo samples
2DMP_AMP_pct <- DMP_AMP %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("A")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

2DMP_FMP_pct <- DMP_FMP %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("M")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

2DMT_AMF_pct <- DMT_AMF %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("A")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

2DMT_FETB_pct <- DMT_FETB %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("F")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

2TDMB_AMP_pct <- TDMB_AMP %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("A")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

2TDMB_FMP_pct <- TDMB_FMP %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("M")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

2TDMT_AMF_pct <- TDMT_AMF %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("A")) <= threshold) * 100) %>%
  ungroup() %>%
   select(gene, PctSamplesPassing)

2TDMT_FETB_pct <- TDMT_FETB %>%
  rowwise() %>%
  mutate(PctSamplesPassing = mean(c_across(contains("F")) <= threshold) * 100) %>%
  ungroup() %>%
  select(gene, PctSamplesPassing)

# check dimensions of in vivo filtered tables
dim(2DMP_AMP_pct)
dim(2DMP_FMP_pct)
dim(2DMT_AMF_pct)
dim(2DMT_FETB_pct)
dim(2TDMB_AMP_pct)
dim(2TDMB_FMP_pct)
dim(2TDMT_AMF_pct)
dim(2TDMT_FETB_pct)


# percentage of genes where more than 80% of samples pass the filter
2more80DMP_AMP <- 2DMP_AMP_pct %>% filter(PctSamplesPassing >80)
dim(2more80DMP_AMP)

2more80DMP_FMP <- 2DMP_FMP_pct %>% filter(PctSamplesPassing >80)
dim(2more80DMP_FMP)

2more80DMT_AMF <- 2DMT_AMF_pct %>% filter(PctSamplesPassing >80)
dim(2more80DMT_AMF)

2more80DMT_FETB <- 2DMT_FETB_pct %>% filter(PctSamplesPassing >80)
dim(2more80DMT_FETB)

2more80TDMB_AMP <- 2TDMB_AMP_pct %>% filter(PctSamplesPassing >80)
dim(2more80TDMB_AMP)

2more80TDMB_FMP <- 2TDMB_FMP_pct %>% filter(PctSamplesPassing >80)
dim(2more80TDMB_FMP)

2more80TDMT_AMF <- 2TDMT_AMF_pct %>% filter(PctSamplesPassing >80)
dim(2more80TDMT_AMF)

2more80TDMT_FETB <- 2TDMT_FETB_pct %>% filter(PctSamplesPassing >80)
dim(2more80TDMT_FETB)

# percentage of genes where more than 80% of samples pass the filter
2more60DMP_AMP <- 2DMP_AMP_pct %>% filter(PctSamplesPassing >60)
dim(2more60DMP_AMP)

2more60DMP_FMP <- 2DMP_FMP_pct %>% filter(PctSamplesPassing >60)
dim(2more60DMP_FMP)

2more60DMT_AMF <- 2DMT_AMF_pct %>% filter(PctSamplesPassing >60)
dim(2more60DMT_AMF)

2more60DMT_FETB <- 2DMT_FETB_pct %>% filter(PctSamplesPassing >60)
dim(2more60DMT_FETB)

2more60TDMB_AMP <- 2TDMB_AMP_pct %>% filter(PctSamplesPassing >60)
dim(2more60TDMB_AMP)

2more60TDMB_FMP <- 2TDMB_FMP_pct %>% filter(PctSamplesPassing >60)
dim(2more60TDMB_FMP)

2more60TDMT_AMF <- 2TDMT_AMF_pct %>% filter(PctSamplesPassing >60)
dim(2more60TDMT_AMF)

2more60TDMT_FETB <- 2TDMT_FETB_pct %>% filter(PctSamplesPassing >60)
dim(2more60TDMT_FETB)

# median percentage of samples passing the filter
2medianDMP_AMP <- 2DMP_AMP_pct %>% summarise(median(PctSamplesPassing))
2medianDMP_AMP

2medianDMP_FMP <- 2DMP_FMP_pct %>% summarise(median(PctSamplesPassing))
2medianDMP_FMP

2medianDMT_AMF <- 2DMT_AMF_pct %>% summarise(median(PctSamplesPassing))
2medianDMT_AMF

2medianDMT_FETB <- 2DMT_FETB_pct %>% summarise(median(PctSamplesPassing))
2medianDMT_FETB

2medianTDMB_AMP <- 2TDMB_AMP_pct %>% summarise(median(PctSamplesPassing))
2medianTDMB_AMP

2medianTDMB_FMP <- 2TDMB_FMP_pct %>% summarise(median(PctSamplesPassing))
2medianTDMB_FMP

2medianTDMT_AMF <- 2TDMT_AMF_pct %>% summarise(median(PctSamplesPassing))
2medianTDMT_AMF

2medianTDMT_FETB <- 2TDMT_FETB_pct %>% summarise(median(PctSamplesPassing))
2medianTDMT_FETB

# density plots individual 
2DMP_AMPdensity <- ggplot(2DMP_AMP_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples Gene (DMP_AMP)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

2DMP_FMPdensity <- ggplot(2DMP_FMP_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (DMP_FMP)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

2DMT_AMFdensity <- ggplot(2DMT_AMF_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (DMT_AMF)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

2DMT_FETBdensity <- ggplot(2DMT_FETB_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (DMT_FETB)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

2TDMB_AMPdensity <- ggplot(2TDMB_AMP_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (TDMB_AMP)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

2TDMB_FMPdensity <- ggplot(2TDMB_FMP_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (TDMB_FMP)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

2TDMT_AMFdensity <- ggplot(2TDMT_AMF_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (TDMT_AMF)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()

2TDMT_FETBdensity <- ggplot(2TDMT_FETB_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "skyblue", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene (TDMT_FETB)",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()



# ecdf plots individual
2DMP_AMPecdf <- ggplot(2DMP_AMP_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (DMP_AMP)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

2DMP_FMPecdf <- ggplot(2DMP_FMP_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (DMP_FMP)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

2DMT_AMFecdf <- ggplot(2DMT_AMF_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (DMT_AMF)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

2DMT_FETBecdf <- ggplot(2DMT_FETB_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (DMT_FETB)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

2TDMB_AMPecdf <- ggplot(2TDMB_AMP_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (TDMB_AMP)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

2TDMB_FMPecdf <- ggplot(2TDMB_FMP_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (TDMB_FMP)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

2TDMT_AMFecdf <- ggplot(2TDMT_AMF_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (TDMT_AMF)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

2TDMT_FETBecdf <- ggplot(2TDMT_FETB_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene (TDMT_FETB)",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()

# combining in vivo filtered tables
2Combined_pct <- bind_rows(DMP_AMP_pct, DMP_FMP_pct, DMT_AMF_pct, DMT_FETB_pct, TDMB_AMP_pct, TDMB_FMP_pct, TDMT_AMF_pct, TDMT_FETB_pct)
dim(2Combined_pct)
head(2Combined_pct)

# percentage of genes where more than 80% of samples pass the filter
2more80 <- 2Combined_pct %>% filter(PctSamplesPassing >80)
2more80

# percentage of genes where more than 60% of samples pass the filter
2more60 <- 2Combined_pct %>% filter(PctSamplesPassing >60)
2more60

# median percentage of samples passing the filter
2medianPCtpassing <- 2Combined_pct %>% summarise(median(PctSamplesPassing))
2medianPCtpassing

# density plot combined in vivo samples
2Combined_densitY <- ggplot(2Combined_pct, aes(x = PctSamplesPassing)) +
  geom_density(fill = "pink", alpha = 0.5) +
  labs(title = "Distribution of % Samples per Gene",
       x = "% of Samples Passing", y = "Density") +
  theme_minimal()
2Combined_density

# combined in vivo tables ecdf
2Combined_ecdf <- ggplot(2Combined_pct, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()
2Combined_ecdf

# combining the in vivo and in vitro filtered tables
Combined1_2 <- bind_rows(1Combined_pct, 2Combined_pct)
dim(Combined1_2)

# total percentage of genes where more than 80% of samples pass the respective filter
more80combi <- Combined1_2 %>% filter(PctSamplesPassing >80)
more80combi

# total percentage of genes where more than 60% of samples pass the respective filter
more60combi <- Combined1_2 %>% filter(PctSamplesPassing >60)
more60combi

# median percenatge of samples passing the respective filter
medianCombi <- Combined1_2 %>% summarise(median(PctSamplesPassing))
medianCombi

# ECDF of combined in vivo and in vitro tables
Combined_ECDFD1_2 <- ggplot(Combined1_2, aes(x = PctSamplesPassing)) +
  stat_ecdf(geom = "step", color = "darkgreen", size = 1) +
  labs(title = "ECDF of % Samples per Gene",
       x = "% of Samples Passing", y = "Cumulative Proportion of Genes") +
  theme_minimal()
Combined_ECDFD1_2
```
[image](/images
