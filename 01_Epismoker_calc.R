####################################################################################

### Calculate the epismoker measure for GS

####################################################################################

## Wave 3

cd /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/00_Preps/

screen

R

W3 = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave3_mvals.rds")
W3[which(is.nan(W3))] <- NA
W3[which(is.infinite(W3))] <- NA

# Transpose so that CpGs are columns
W3 <- t(W3)

# Convert to beta 
meth <- W3
m2beta <- function(meth) { 
  beta <- 2^meth/(2^meth + 1)
  return(beta)
}
W3 <- m2beta(meth) # convert m-value to beta-value 

# Impute missing
library(imputeTS)
W3 <- na_mean(W3)


W3 <- t(W3) # cpgs as rows needed 

# Now apply EpiSmoker to calculate scores 
library(EpiSmokEr)
result <- epismoker(dataset=W3, method = "SSc")

# Save out data scores for smoking
saveRDS(result, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave3_epismoker_updated_cpgs_as_rows.rds")


## Wave 4

cd /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/00_Preps/

screen

R

W4 = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/wave4/w4-mvals.rds")

# Transpose so that CpGs are columns
W4 <- t(W4)

# Convert to beta 
meth <- W4
m2beta <- function(meth) { 
  beta <- 2^meth/(2^meth + 1)
  return(beta)
}
W4 <- m2beta(meth) # convert m-value to beta-value 

# Impute missing
library(imputeTS)
W4 <- na_mean(W4)

W4 <- t(W4) # cpgs as rows needed 

# Now apply EpiSmoker to calculate scores 
library(EpiSmokEr)
result <- epismoker(dataset=W4, method = "SSc")

# Save out data scores for smoking
saveRDS(result, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave4_epismoker_updated_cpgs_as_rows.rds")


## Wave 1

cd /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/00_Preps/

screen

R

meth = readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/norm_mvals_5087.rds")

# Transpose so that CpGs are columns
meth <- t(meth)

# Convert to beta 
m2beta <- function(meth) { 
  beta <- 2^meth/(2^meth + 1)
  return(beta)
}
meth <- m2beta(meth) # convert m-value to beta-value 

# Impute missing
library(imputeTS)
meth <- na_mean(meth)

meth <- t(meth) # cpgs as rows needed 

# Now apply EpiSmoker to calculate scores 
library(EpiSmokEr)
result <- epismoker(dataset=meth, method = "SSc")

# Save out data scores for smoking
saveRDS(result, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave1_epismoker_updated_cpgs_as_rows.rds")


###################################################

### Correlate EpiSmokR with pack years in GS as a sense check 

# Read in epismoekr 
w1 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave1_epismoker_updated_cpgs_as_rows.rds")
w3 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave3_epismoker_updated_cpgs_as_rows.rds")
w4 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave4_epismoker_updated_cpgs_as_rows.rds")

bind <- rbind(w1, w3)
bind <- rbind(bind, w4) # 18779 individuals with DNAm epismoker calculated 
bind$Sample_Sentrix_ID <- row.names(bind)

# Read in pack years info 

packyears <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/updated_smoking_jan_2019/pack_years.csv")
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
library(tidyverse)
packyears <- left_join(packyears, target, by = "Sample_Name")


merge <- merge(bind, packyears, by = "Sample_Sentrix_ID") # 18413

merge <- merge[complete.cases(merge$pack_years),]

cor.test(merge$pack_years, merge$smokingScore)


merge_w1 <- merge[which(merge$Set == "wave1"),]

cor.test(merge_w1$pack_years, merge_w1$smokingScore)



####################################################################################

### Calculate the epismoker measure for LBC1936

####################################################################################

source("http://bioconductor.org/biocLite.R")
install.packages("devtools") # if you don't have the package, run install.packages("devtools")
library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")
install.packages("IlluminaHumanMethylation450kmanifest")

library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)

suppressPackageStartupMessages({
library(EpiSmokEr)  
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(htmlTable)
library(rmarkdown)
})

Read in target file with matched IDs for methylation data 
load("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/Beta_3525_norm_bgcorrect_0.001BetaThreshold_probefilter.RObject")
result2 <- epismoker(dataset=dat, method = "SSc")

Save out data scores for smoking
write.rds("/Cluster_Filespace/Marioni_Group/Danni/lbc_epismoker.rds") 

