#############################################################################################################

### Protein and DNAm preps - GDF15 and Nt-pro-BNP

#############################################################################################################

# Protein measures are first prepared for running in elent models
# Next, DNAm is prepped in the relevant individuals for train/test
# These inputs then feed into the elnet scripts for predictor score generation
# DNAm is prepped for 450k array and EPIC array 

#############################################################################################################

### Prep GDF15 and Nt-pro-BNP test sets

#############################################################################################################

screen

R

library(tidyverse)
library(readxl)

# Read in variable data for proteins
d1 <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/file_input_030821/GS20K_GDF15_NT_proBNP.PHE")
names(d1)[2] <- "Sample_Name"

# Read in GS 10k meth target file and join into phenos
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
d1 <- left_join(d1, target, by = "Sample_Name")

# Get complete pheno info for each biomarker separately and joint
table(is.na(d1$gdf15)) # 18414
table(is.na(d1$nt.probnp)) # 17863 

GDF <- d1 %>% filter(gdf15 != "NA")
BNP <- d1 %>% filter(nt.probnp != "NA")

### PREP EACH BIOMARKER SEPARATELY

## GDF 

# Assign biomarker name to these parameters 
prot <- GDF 
name <- "GDF"
biomarker <- "gdf15"

# Join target and restrict to those with DNAm data available - then save a copy 
d1 <- prot %>% filter(Batch != "NA") # 17489

# # Split into those in waves 1 and 3 of GS to see how many we have that should match 
W1 <- d1 %>% filter(Set == "wave1")
W3 <- d1 %>% filter(Set == "wave3")
W4 <- d1 %>% filter(Set == "wave4")

# > dim(W4)
# [1] 8369   33
# > dim(W3)
# [1] 4315   33
# > dim(W1)
# [1] 4805   33

# Read in genetic pQTL information table and sentinel pQTL data file from GWAS lookup step
table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/pQTL_extraction_index_table.csv")
BNP_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_formatted/all_BNP.csv", check.names = F)
GDF_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/GDF_formatted/GDF15_chr_19.csv", check.names = F)

table <- table %>% filter(Protein == "GDF15")

# We need to join the 3 extracted SNPs into the d1 file for residualisation
names(GDF_ex)[1] <- "Sample_Name"
d1 <- left_join(d1, GDF_ex, by = "Sample_Name")

# Residualise phenotypes in separate train/test sets for elnets
W1 <- left_join(W1, GDF_ex, by = "Sample_Name")
W3 <- left_join(W3, GDF_ex, by = "Sample_Name")
W4 <- left_join(W4, GDF_ex, by = "Sample_Name")

# Train and test sets 
W1 <- rbind(W1,W3) # 9120 GDF

# Exclusions - remove relatedness from training set 
GDF_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/GDF_exclusion_training_W1W3.csv")
BNP_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/BNP_exclusion_training_W1W3.csv")

# Exclusions - keep these for test set 
GDF_keep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/unrelated_W4_GDF.csv")
BNP_keep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/unrelated_W4_BNP.csv")

W1 <- W1[-which(W1$Sample_Name %in% GDF_ex$Sample_Name),]
W4 <- W4[which(W4$Sample_Name %in% GDF_keep$Sample_Name),]

# > dim(W1)
# [1] 8207   36
# > dim(W4)
# [1] 2954   36

# New updated: 
# > dim(W1)
# [1] 8362   36
# > dim(W4)
# [1] 2622   36


## Rank-Inverse Based Normaliation 
library(bestNormalize)
for(i in colnames(W1)[c(4,6)]){ 
  W1[,i]<- orderNorm(W1[,i])$x.t
}

## Rank-Inverse Based Normaliation 
library(bestNormalize)
for(i in colnames(W4)[c(4,6)]){ 
  W4[,i]<- orderNorm(W4[,i])$x.t
}

# Residualise in train and test sets and scale
list <- c(4,6)
for(i in list) {
    W1[,i] <- scale(resid(lm(W1[,i] ~ age.x + factor(sex.x) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
         + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data=W1, na.action=na.exclude)))
}

# list <- c(4,6)
# for(i in list) {
#     W4[,i] <- scale(resid(lm(W4[,i] ~ age.x + factor(sex.x) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
#          + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data=W4, na.action=na.exclude)))
# }

# Write out train and test proteins without pQTL adjustment, but with asjustment for residualisation and scaling/transforming
W4 <- W4[c(29,4)]
write.csv(W4, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W4_GDF15.csv", row.names = F)

W1 <- W1[c(29,4)]
write.csv(W1, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W1_GDF15.csv", row.names = F)


###################

## BNP

# Assign biomarker name to these parameters 
prot <- BNP
name <- "BNP"
biomarker <- "nt.probnp"

# Join target and restrict to those with DNAm data available - then save a copy 
d1 <- prot %>% filter(Batch != "NA") # 16963
write.csv(d1, paste0("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/01_elnets_20k/00_Preps/", name,"_16963_with_DNAm_phen_reference_file.csv"), row.names = F)

# # Split into those in waves 1 and 3 of GS to see how many we have that should match 
W1 <- d1 %>% filter(Set == "wave1")
W3 <- d1 %>% filter(Set == "wave3")
W4 <- d1 %>% filter(Set == "wave4")

# > dim(W1)
# [1] 4662   33
# > dim(W3)
# [1] 4178   33
# > dim(W4)
# [1] 8123   33

# Read in genetic pQTL information table and sentinel pQTL data file from GWAS lookup step
table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/pQTL_extraction_index_table.csv")
BNP_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_formatted/all_BNP.csv", check.names = F)
GDF_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/GDF_formatted/GDF15_chr_19.csv", check.names = F)

table <- table %>% filter(Protein == "NtproBNP")

# We need to join the 10 extracted SNPs into the d1 file for residualisation
names(BNP_ex)[1] <- "Sample_Name"
d1 <- left_join(d1, BNP_ex, by = "Sample_Name")

# Residualise phenotypes in separate train/test sets for elnets
W1 <- left_join(W1, BNP_ex, by = "Sample_Name")
W3 <- left_join(W3, BNP_ex, by = "Sample_Name")
W4 <- left_join(W4, BNP_ex, by = "Sample_Name")

# Train and test sets 
W1 <- rbind(W1,W3) # 8840 BNP

# Exclusions - remove relatedness from training set 
GDF_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/GDF_exclusion_training_W1W3.csv")
BNP_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/BNP_exclusion_training_W1W3.csv")

# Exclusions - keep these for test set 
GDF_keep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/unrelated_W4_GDF.csv")
BNP_keep <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/unrelated_W4_BNP.csv")

W1 <- W1[-which(W1$Sample_Name %in% BNP_ex$Sample_Name),]
W4 <- W4[which(W4$Sample_Name %in% BNP_keep$Sample_Name),]


# > dim(W1)
# [1] 8144   43
# > dim(W4)
# [1] 2476   43


## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes
library(bestNormalize)
for(i in colnames(W1)[c(4,6)]){ 
  W1[,i]<- orderNorm(W1[,i])$x.t
}

## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes
library(bestNormalize)
for(i in colnames(W4)[c(4,6)]){ 
  W4[,i]<- orderNorm(W4[,i])$x.t
}

# Residualise in train and test and scale
list <- c(4,6)
for(i in list) {
    W1[,i] <- scale(resid(lm(W1[,i] ~ age.x + factor(sex.x) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
         + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data=W1, na.action=na.exclude)))
}


# list <- c(4,6)
# for(i in list) {
#     W4[,i] <- scale(resid(lm(W4[,i] ~ age.x + factor(sex.x) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
#          + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data=W4, na.action=na.exclude)))
# }


# Write out train and test proteins without pQTL adjustment, but with asjustment for residualisation and scaling/transforming
W4 <- W4[c(29,6)]
write.csv(W4, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W4_BNP.csv", row.names = F)

W1 <- W1[c(29,6)]
write.csv(W1, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W1_BNP.csv", row.names = F)


#############################################################################################################

### Prep GDF15 and Nt-pro-BNP levels - for elnets in 20k full GS 

#############################################################################################################

### RESIDUALISE BY WAVES AND COMBINE TO ELNET TRAIN/TEST FILES
# Regress onto pQTLs for each biomarker

screen

R

library(tidyverse)
library(readxl)

# Read in variable data for proteins
d1 <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/file_input_030821/GS20K_GDF15_NT_proBNP.PHE")
names(d1)[2] <- "Sample_Name"

# Read in GS 10k meth target file and join into phenos
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
d1 <- left_join(d1, target, by = "Sample_Name")

# Get complete pheno info for each biomarker separately and joint
table(is.na(d1$gdf15)) # 18414
table(is.na(d1$nt.probnp)) # 17863 

GDF <- d1 %>% filter(gdf15 != "NA")
BNP <- d1 %>% filter(nt.probnp != "NA")

### PREP EACH BIOMARKER SEPARATELY

## GDF 

# Assign biomarker name to these parameters 
prot <- GDF 
name <- "GDF"
biomarker <- "gdf15"

# Join target and restrict to those with DNAm data available - then save a copy 
d1 <- prot %>% filter(Batch != "NA") # 17489

# # Split into those in waves 1 and 3 of GS to see how many we have that should match 
W1 <- d1 %>% filter(Set == "wave1")
W3 <- d1 %>% filter(Set == "wave3")
W4 <- d1 %>% filter(Set == "wave4")

# > dim(W4)
# [1] 8369   33
# > dim(W3)
# [1] 4315   33
# > dim(W1)
# [1] 4805   33

# Read in genetic pQTL information table and sentinel pQTL data file from GWAS lookup step
table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/pQTL_extraction_index_table.csv")
BNP_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_formatted/all_BNP.csv", check.names = F)
GDF_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/GDF_formatted/GDF15_chr_19.csv", check.names = F)

table <- table %>% filter(Protein == "GDF15")

# We need to join the 3 extracted SNPs into the d1 file for residualisation
names(GDF_ex)[1] <- "Sample_Name"
d1 <- left_join(d1, GDF_ex, by = "Sample_Name")

# Residualise phenotypes in separate train/test sets for elnets
W1 <- left_join(W1, GDF_ex, by = "Sample_Name")
W3 <- left_join(W3, GDF_ex, by = "Sample_Name")
W4 <- left_join(W4, GDF_ex, by = "Sample_Name")

# Train 
W1 <- rbind(W1,W3) 
W1 <- rbind(W1,W4) # 17489

## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes
library(bestNormalize)
for(i in colnames(W1)[c(4,6)]){ 
  W1[,i]<- orderNorm(W1[,i])$x.t
}


# Residualise in train 
list <- c(4,6)
for(i in list) {
    W1[,i] <- scale(resid(lm(W1[,i] ~ age.x + factor(sex.x) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
         + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data=W1, na.action=na.exclude)))
}


# Write out train and test proteins with pQTL adjustment and with asjustment for residualisation and scaling/transforming
W1 <- W1[c(29,4)]
write.csv(W1, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/ALL_GDF15.csv", row.names = F)


###################

## BNP

# Assign biomarker name to these parameters 
prot <- BNP
name <- "BNP"
biomarker <- "nt.probnp"

# Join target and restrict to those with DNAm data available - then save a copy 
d1 <- prot %>% filter(Batch != "NA") # 16963
write.csv(d1, paste0("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/01_elnets_20k/00_Preps/", name,"_16963_with_DNAm_phen_reference_file.csv"), row.names = F)

# # Split into those in waves 1 and 3 of GS to see how many we have that should match 
W1 <- d1 %>% filter(Set == "wave1")
W3 <- d1 %>% filter(Set == "wave3")
W4 <- d1 %>% filter(Set == "wave4")

# > dim(W1)
# [1] 4662   33
# > dim(W3)
# [1] 4178   33
# > dim(W4)
# [1] 8123   33

# Read in genetic pQTL information table and sentinel pQTL data file from GWAS lookup step
table <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/pQTL_extraction_index_table.csv")
BNP_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/BNP_formatted/all_BNP.csv", check.names = F)
GDF_ex <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/pQTL_extraction/GDF_formatted/GDF15_chr_19.csv", check.names = F)

table <- table %>% filter(Protein == "NtproBNP")

# We need to join the 10 extracted SNPs into the d1 file for residualisation
names(BNP_ex)[1] <- "Sample_Name"
d1 <- left_join(d1, BNP_ex, by = "Sample_Name")

# Residualise phenotypes in separate train/test sets for elnets
W1 <- left_join(W1, BNP_ex, by = "Sample_Name")
W3 <- left_join(W3, BNP_ex, by = "Sample_Name")
W4 <- left_join(W4, BNP_ex, by = "Sample_Name")

# Train and test sets 
W1 <- rbind(W1,W3) 
W1 <- rbind(W1,W4) # 16963

## Rank-Inverse Based Normaliation on Pre-corrected Protein Phenotypes
library(bestNormalize)
for(i in colnames(W1)[c(4,6)]){ 
  W1[,i]<- orderNorm(W1[,i])$x.t
}

# Residualise in train and test and scale
list <- c(4,6)
for(i in list) {
    W1[,i] <- scale(resid(lm(W1[,i] ~ age.x + factor(sex.x) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
         + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data=W1, na.action=na.exclude)))
}


W1 <- W1[c(29,6)]
write.csv(W1, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/ALL_BNP.csv", row.names = F)



#############################################################################################################

### PREP DNAm 450k data separately and save out files for use in elnets

#############################################################################################################

cd /Local_Scratch/Danni/GDFBNP/00_Preps/

screen

R

### Prep each wave for 450k array in beta vals

## Wave 1
W1 = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/norm_mvals_5087.rds")

# Transpose so that CpGs are columns
W1 <- t(W1)

# # Subset to 450k
# anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
# common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]
# W1 <- W1[,colnames(W1) %in% rownames(common_anno)] 

# Convert to beta 
m2beta <- function(W1) { 
  beta <- 2^W1/(2^W1 + 1)
  return(beta)
}
W1 <- m2beta(W1) # convert m-value to beta-value 

# Impute missing
library(imputeTS)
W1 <- na_mean(W1)

# Remove cross-hyb, snp-associated and XY cpgs
library(tidyverse)
snps <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/snp_probes.txt", sep='\t', header=T)
snps <- snps[which(snps$EUR_AF >= 0.05), "IlmnID"] %>% as.character
ch1 <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cph.txt", sep='\t', header=F)
ch1 <- as.character(ch1[,1])
ch2 <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cpg.txt", sep='\t', header=F)
ch2 <- as.character(ch2[,1])
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_UCSC/annotations_UCSC_object_df_from_daniel.rds")
ychr <- anno[anno$chr=="chrY", "Name"]
exclude <- c(snps, ch1, ch2, ychr) %>% unique # 54,592
W1 <- W1[,-which(colnames(W1) %in% exclude)]
dim(W1)


## Wave 3
W3 = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3_mvals.rds")
W3[which(is.nan(W3))] <- NA
W3[which(is.infinite(W3))] <- NA

# Transpose so that CpGs are columns
W3 <- t(W3)

# # Subset to 450k
# anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
# common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]
# W3 <- W3[,colnames(W3) %in% rownames(common_anno)] #  4450 398624

# Convert to beta 
m2beta <- function(W3) { 
  beta <- 2^W3/(2^W3 + 1)
  return(beta)
}
W3 <- m2beta(W3) # convert m-value to beta-value 

# Impute missing
library(imputeTS)
W3 <- na_mean(W3)

# Remove cross-hyb, snp-associated and XY cpgs
library(tidyverse)
snps <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/snp_probes.txt", sep='\t', header=T)
snps <- snps[which(snps$EUR_AF >= 0.05), "IlmnID"] %>% as.character
ch1 <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cph.txt", sep='\t', header=F)
ch1 <- as.character(ch1[,1])
ch2 <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cpg.txt", sep='\t', header=F)
ch2 <- as.character(ch2[,1])
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_UCSC/annotations_UCSC_object_df_from_daniel.rds")
ychr <- anno[anno$chr=="chrY", "Name"]
exclude <- c(snps, ch1, ch2, ychr) %>% unique # 54,592
W3 <- W3[,-which(colnames(W3) %in% exclude)]
dim(W3) # 4450 393654


## Wave 4
W4 = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/w4-mvals.rds")

# Transpose so that CpGs are columns
W4 <- t(W4)

# Subset to individuals
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/w4-samplesheet_v3.rds")
W4 <- W4[which(rownames(W4) %in% rownames(target)),] # 8877 854642

# # Load EPIC array file and subset to probes common to 450k and EPIC array - across test and train sets 
# anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
# common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]
# W4 <- W4[,colnames(W4) %in% rownames(common_anno)] #  8877 445235

# Convert to beta 
m2beta <- function(W4) { 
  beta <- 2^W4/(2^W4 + 1)
  return(beta)
}
W4 <- m2beta(W4) # convert m-value to beta-value 

# Impute missing
library(imputeTS)
W4 <- na_mean(W4)

# Remove cross-hyb, snp-associated and XY cpgs
library(tidyverse)
snps <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/snp_probes.txt", sep='\t', header=T)
snps <- snps[which(snps$EUR_AF >= 0.05), "IlmnID"] %>% as.character
ch1 <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cph.txt", sep='\t', header=F)
ch1 <- as.character(ch1[,1])
ch2 <- read.table("/Cluster_Filespace/Marioni_Group/mQTL_and_Illumina_annotation_files/epic_probes_to_exclude/nonspecific_cpg.txt", sep='\t', header=F)
ch2 <- as.character(ch2[,1])
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/Stradl_markers/00_Revisions_updates/Annotations_UCSC/annotations_UCSC_object_df_from_daniel.rds")
ychr <- anno[anno$chr=="chrY", "Name"]
exclude <- c(snps, ch1, ch2, ychr) %>% unique # 54,592
W4 <- W4[,-which(colnames(W4) %in% exclude)]
dim(W4) # 8877 411979

###################################################################

## Ensure probes common to all waves are retained

# > dim(W1)
# [1]   5087 807785
# > dim(W3)
# [1]   4450 765695
# > dim(W4)
# [1]   8877 802624

# Get common probes across all waves of train/test and subset 
W3 <- W3[,which(colnames(W3) %in% colnames(W1))]
W3 <- W3[,which(colnames(W3) %in% colnames(W4))]

W1 <- W1[,which(colnames(W1) %in% colnames(W3))]
W1 <- W1[,which(colnames(W1) %in% colnames(W4))]

W4 <- W4[,which(colnames(W4) %in% colnames(W1))]
W4 <- W4[,which(colnames(W4) %in% colnames(W3))]

# > dim(W1)
# [1]   5087 760838
# > dim(W3)
# [1]   4450 760838
# > dim(W4)
# [1]   8877 760838

## Make sure order of CpG sites match in files for joining  
id = colnames(W3)
W1 <- W1[,match(id, colnames(W1))]
table(colnames(W3) == colnames(W1))

id = colnames(W3)
W4 <- W4[,match(id, colnames(W4))]
table(colnames(W3) == colnames(W4))
table(colnames(W1) == colnames(W4))

# Save out train and test sets 
library(data.table)
saveRDS(W1, "/Local_Scratch/Danni/GDFBNP/00_Preps/W1_EPIC.rds")

saveRDS(W3, "/Local_Scratch/Danni/GDFBNP/00_Preps/W3_EPIC.rds")

saveRDS(W4, "/Local_Scratch/Danni/GDFBNP/00_Preps/W4_EPIC.rds")

##############################################################################

### Subset DNAm to create x files based on y files (450k)

cd /Local_Scratch/Danni/GDFBNP/00_Preps/

screen

R

# Read in wave DNAm prepped
W1 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/W1_EPIC.rds")
W3 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/W3_EPIC.rds")
W4 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/W4_EPIC.rds")

# > dim(W1)
# [1]   5087 760838
# > dim(W3)
# [1]   4450 760838
# > dim(W4)
# [1]   8877 760838

# y files (taken from preps above)
BNP_train <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W1_BNP.csv")
BNP_test <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W4_BNP.csv")

GDF_train <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W1_GDF15.csv")
GDF_test <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W4_GDF15.csv")

# > dim(BNP_train)
# [1] 8002    2
# > dim(BNP_test)
# [1] 2808    2
# > dim(GDF_train)
# [1] 8207    2
# > dim(GDF_test)
# [1] 2954    2

## GDF

# GDF w1
GDF_W1 <- W1[which(rownames(W1) %in% GDF_train$Sample_Sentrix_ID),]
# GDF w3
GDF_W3 <- W3[which(rownames(W3) %in% GDF_train$Sample_Sentrix_ID),]
# GDF w4
GDF_W4 <- W4[which(rownames(W4) %in% GDF_test$Sample_Sentrix_ID),] 

W1 = NULL
W3 = NULL
W4 = NULL

# Scale CpGs within wave
names <- rownames(GDF_W1)
GDF_W1 <- apply(GDF_W1, 2, scale)
rownames(GDF_W1) <- names

names <- rownames(GDF_W3)
GDF_W3 <- apply(GDF_W3, 2, scale)
rownames(GDF_W3) <- names

# Combine W1 W3 for training 
GDF_W1W3 <- rbind(GDF_W1, GDF_W3)

GDF_W1 = NULL
GDF_W3 = NULL

names <- rownames(GDF_W4)
GDF_W4 <- apply(GDF_W4, 2, scale)
rownames(GDF_W4) <- names



# Ensure sample order matches between x and y files for training and testing
GDF_W1W3 <- GDF_W1W3[match(GDF_train$Sample_Sentrix_ID, rownames(GDF_W1W3)),]
identical(GDF_train$Sample_Sentrix_ID, rownames(GDF_W1W3))

GDF_W4 <- GDF_W4[match(GDF_test$Sample_Sentrix_ID, rownames(GDF_W4)),]
identical(GDF_test$Sample_Sentrix_ID, rownames(GDF_W4))

# Save out x files EPIC
library(data.table)
saveRDS(GDF_W1W3, "/Local_Scratch/Danni/GDFBNP/00_Preps/GDF_train_x_subset_EPIC.rds")
saveRDS(GDF_W4, "/Local_Scratch/Danni/GDFBNP/00_Preps/GDF_test_x_subset_EPIC.rds")

# Subset to 450k array versions
# # Load EPIC array file and subset to probes common to 450k and EPIC array - across test and train sets 
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

GDF_W4 <- GDF_W4[,colnames(GDF_W4) %in% rownames(common_anno)] # 2622 390461
GDF_W1W3 <- GDF_W1W3[,colnames(GDF_W1W3) %in% rownames(common_anno)] # 8362 390461

# Save out x files 450k
saveRDS(GDF_W1W3, "/Local_Scratch/Danni/GDFBNP/00_Preps/GDF_train_x_subset.rds")
saveRDS(GDF_W4, "/Local_Scratch/Danni/GDFBNP/00_Preps/GDF_test_x_subset.rds")
GDF_W1W3 = NULL
GDF_W4 = NULL


########################################################################################

## BNP

cd /Local_Scratch/Danni/GDFBNP/00_Preps/

screen

R

# Read in wave DNAm prepped
W1 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/W1_EPIC.rds")
W3 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/W3_EPIC.rds")
W4 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/W4_EPIC.rds")

# > dim(W1)
# [1]   5087 760838
# > dim(W3)
# [1]   4450 760838
# > dim(W4)
# [1]   8877 760838

# y files (taken from preps above)
BNP_train <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W1_BNP.csv")
BNP_test <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W4_BNP.csv")

GDF_train <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W1_GDF15.csv")
GDF_test <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W4_GDF15.csv")


## BNP

# BNP w1
BNP_W1 <- W1[which(rownames(W1) %in% BNP_train$Sample_Sentrix_ID),]
# BNP w3
BNP_W3 <- W3[which(rownames(W3) %in% BNP_train$Sample_Sentrix_ID),]
# BNP w4
BNP_W4 <- W4[which(rownames(W4) %in% BNP_test$Sample_Sentrix_ID),] 

W1 = NULL
W3 = NULL
W4 = NULL


# Scale CpGs within wave
names <- rownames(BNP_W1)
BNP_W1 <- apply(BNP_W1, 2, scale)
rownames(BNP_W1) <- names

names <- rownames(BNP_W3)
BNP_W3 <- apply(BNP_W3, 2, scale)
rownames(BNP_W3) <- names

# Combine W1 W3 for training 
BNP_W1W3 <- rbind(BNP_W1, BNP_W3)

BNP_W1 = NULL
BNP_W3 = NULL

names <- rownames(BNP_W4)
BNP_W4 <- apply(BNP_W4, 2, scale)
rownames(BNP_W4) <- names


# Ensure sample order matches between x and y files for training and testing
BNP_W1W3 <- BNP_W1W3[match(BNP_train$Sample_Sentrix_ID, rownames(BNP_W1W3)),]
identical(BNP_train$Sample_Sentrix_ID, rownames(BNP_W1W3))

BNP_W4 <- BNP_W4[match(BNP_test$Sample_Sentrix_ID, rownames(BNP_W4)),]
identical(BNP_test$Sample_Sentrix_ID, rownames(BNP_W4))

# Save out train and test sets 
library(data.table)
saveRDS(BNP_W1W3, "/Local_Scratch/Danni/GDFBNP/00_Preps/BNP_train_x_subset_EPIC.rds")
saveRDS(BNP_W4, "/Local_Scratch/Danni/GDFBNP/00_Preps/BNP_test_x_subset_EPIC.rds")

# Subset to 450k array versions
# # Load EPIC array file and subset to probes common to 450k and EPIC array - across test and train sets 
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

BNP_W4 <- BNP_W4[,colnames(BNP_W4) %in% rownames(common_anno)] #  390461
BNP_W1W3 <- BNP_W1W3[,colnames(BNP_W1W3) %in% rownames(common_anno)] #  390461

# Save out x files 450k
saveRDS(BNP_W1W3, "/Local_Scratch/Danni/GDFBNP/00_Preps/BNP_train_x_subset.rds")
saveRDS(BNP_W4, "/Local_Scratch/Danni/GDFBNP/00_Preps/BNP_test_x_subset.rds")

BNP_W1W3 = NULL
BNP_W4 = NULL

## Now x and y files are ready to feed into elnets for within GS train/test




###################################################################################

### PREP TRAINING SET FOR THE 20K EPISCORE 

###################################################################################

cd /Local_Scratch/Danni/GDFBNP/00_Preps/

screen

R

# Read in wave DNAm prepped
W1 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/W1_EPIC.rds")
W3 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/W3_EPIC.rds")
W4 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/W4_EPIC.rds")

# Scale CpGs within wave
names <- rownames(W1)
W1 <- apply(W1, 2, scale)
rownames(W1) <- names

names <- rownames(W3)
W3 <- apply(W3, 2, scale)
rownames(W3) <- names

names <- rownames(W4)
W4 <- apply(W4, 2, scale)
rownames(W4) <- names

# Combine W1 W3 for training 
W1 <- rbind(W1,W3)
W3 = NULL

W1 <- rbind(W1,W4)
W4 = NULL

# Save out EPIC and 450k versions that have been scaled in separate waves before combining
library(data.table)
saveRDS(W1, "/Local_Scratch/Danni/GDFBNP/00_Preps/20k_EPIC.rds") #  18414 and 760838

# Subset to 450k array versions
# # Load EPIC array file and subset to probes common to 450k and EPIC array
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Daniel/EPIC_AnnotationObject_df.rds")
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

W1 <- W1[,colnames(W1) %in% rownames(common_anno)] # 18414 390461

# Save out x files 450k
saveRDS(W1, "/Local_Scratch/Danni/GDFBNP/00_Preps/20k_450K.rds")

W1 = NULL

## Now x and y files are ready to feed into elnets for within GS train/test (must be subset to exact values in train script and matched)

