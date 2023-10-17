###############################################################################################################

### PREP MWAS FILES

###############################################################################################################

## Open screen and start session in R
# cd /Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/
# screen
# R

## Set working directory 
setwd("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/") 

## Load requisite libraries
library(data.table)
library(limma)
library(lumi)
library(readxl)
library(tidyverse)

## Create function for mean imputation
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

#########################################################
######## STEP 1 - PREPARATION OF PHENOTYPE FILE ########
#########################################################

# Read in new variable data from Riccardo (taken from GS folder on datastore)
d1 <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/file_input_030821/GS20K_GDF15_NT_proBNP.PHE")
names(d1)[2] <- "Sample_Name"

# Get complete pheno info 
table(is.na(d1$nt.probnp_rnk)) # 17863
table(is.na(d1$gdf15_rnk)) # 18414
GDF15 <- d1[complete.cases(d1$gdf15_rnk),]
NT <- d1[complete.cases(d1$nt.probnp),] 

# Read in 20k DNAm targets and add info
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
length(complete.cases(target$Sample_Sentrix_ID)) # 18413
GDF15 <- left_join(GDF15, target, by = "Sample_Name")
NT <- left_join(NT, target, by = "Sample_Name")

# Isolate those with each protein + DNAm
GDF15 <- GDF15[complete.cases(GDF15$Sample_Sentrix_ID),] # 17489
NT <- NT[complete.cases(NT$Sample_Sentrix_ID),] # 16963

# Make reference files with FID and IID for each protein
ids1 <- data.frame(FID = GDF15$FID, IID = GDF15$Sample_Name)
ids2 <- data.frame(FID = NT$FID, IID = NT$Sample_Name)
write.table(ids1,"/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/gdf_ids.list",row.names=F,quote=F,sep="\t")
write.table(ids2,"/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/bnp_ids.list",row.names=F,quote=F,sep="\t")

# Residualise GDF (rank transformed protein)
list <- c(3)
for(i in list) {
    GDF15[,i] <- scale(resid(lm(GDF15[,i] ~ age.x + factor(sex.x) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
         + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data=GDF15, na.action=na.exclude)))
}

# Residualise BNP (rank transformed protein)
list <- c(5)
for(i in list) {
    NT[,i] <- scale(resid(lm(NT[,i] ~ age.x + factor(sex.x) + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10
         + PC11 + PC12 + PC13 + PC14 + PC15 + PC16 + PC17 + PC18 + PC19 + PC20, data=NT, na.action=na.exclude)))
}

# Index ids
id1 = GDF15$Sample_Sentrix_ID
id2 = NT$Sample_Sentrix_ID

# Write out phenotypes
write.table(x = t(as.matrix(as.numeric(GDF15[,3]))),file="/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/gdf_resid.csvphen",quote = F, sep = ",", row.names = F, col.names = F)
write.table(x = t(as.matrix(as.numeric(NT[,5]))),file="/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/bnp_resid.csvphen",quote = F, sep = ",", row.names = F, col.names = F)

# Save protein files to subset against 
write.csv(GDF15, '/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/GDF15_protein_file.csv', row.names = F)
write.csv(NT, '/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/BNP_protein_file.csv', row.names = F)

##########################################################
######## STEP 2 - PREPARATION OF METHYLATION FILE GDF ####
##########################################################

## Open screen and start session in R
# cd /Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/
# screen
# R

## Load requisite libraries
library(data.table)
library(limma)
library(lumi)
library(readxl)
library(tidyverse)

## Create function for mean imputation
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

setwd("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/") 

# Load phenotype file for GDF and process DNAm
phenos <- read.csv("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/GDF15_protein_file.csv")

# Read in methylation basenames  
samps_ref=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds") # 18413
samps_ref <- samps_ref[which(samps_ref$Sample_Name %in% phenos$Sample_Name),] # 17489

# Read in methylation file # by chromosome 
for(i in 1:22){ 
  meth=readRDS(paste0("GS20k_chr", i, "_mvals.rds"))
  
  # Subset to those with phenotype in question 
  meth=meth[,which(colnames(meth) %in% phenos$Sample_Sentrix_ID)] #80545 17464
  
  # Subset to probes passing QC 
  probes=read.table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)
  meth=meth[which(row.names(meth) %in% probes$V1),] # 73108 17464
  
  # # Subset to 450k only 
  # meth=meth[which(row.names(meth)%in%anno$Name),]
  
  # Match order of IDs in phenotype and methylation file 
  ids=phenos$Sample_Sentrix_ID
  meth=meth[,match(ids,colnames(meth))]
  
  # Check order of IDs match between phenotype and methylation files 
  print(table(colnames(meth)==phenos$Sample_Sentrix_ID))
  
  # Convert to beta values 
  meth=m2beta(meth)
  
  # Mean impute - cannot have NAs in final file 
  meth <- apply(meth,1,meanimpute)
  
  # Transpose back original format - as apply changes the format
  meth=t(meth)
  
  # Prepare covariate matrix for regressions 
  # Match order of IDs with other files 
  samps <- samps_ref
  samps=samps[match(ids,samps$Sample_Sentrix_ID),] # 17464
  
  # Check order of IDs match with other files 
  print(table(samps$Sample_Sentrix_ID==phenos$Sample_Sentrix_ID))
  print(table(samps$Sample_Sentrix_ID==colnames(meth)))
  
  ## Regression step - residualise for age, sex and batch 
  design.resid <- model.matrix(~sex + age + Batch, data=samps)
  fit.resid <- limma::lmFit(meth, design.resid)
  gc()
  meth <- limma::residuals.MArrayLM(fit.resid, meth)
  meth <- meth[!is.infinite(rowSums(meth)), ]
  rm(fit.resid)
  gc()
  
  ## Scale 
  meth=t(apply(meth,1,scale))
  
  ## Write out CpGs 
  cpgs=as.data.frame(row.names(meth))
  names(cpgs)[1]="CpG"
  fwrite(cpgs, paste0("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/cpgs_chr_index_gdf/GS20k_chr", i, "_subset_cpgs.txt"),row.names=F)
  
  # Save out residualised file 
  fwrite(meth, paste0("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/meth_resid_gdf/GS20k_chr", i, "_resid.txt"),row.names=F)
  
  ## Remove methylation object and clean up environment 
  rm(meth)
  gc()
  print(i)
} 

## Write out ID order for later
saveRDS(ids, "/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/gdf_meth_ids.rds")
 

##########################################################
######## STEP 2 - PREPARATION OF METHYLATION FILE BNP ####
##########################################################

## Open screen and start session in R
# cd /Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/
# screen
# R

## Load requisite libraries
library(data.table)
library(limma)
library(lumi)
library(readxl)
library(tidyverse)

## Create function for mean imputation
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

setwd("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/Chromosomes/") 

# Load phenotype file for GDF and process DNAm
phenos <- read.csv("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/BNP_protein_file.csv")

# Read in methylation basenames  
samps_ref=readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds") # 18413
samps_ref <- samps_ref[which(samps_ref$Sample_Name %in% phenos$Sample_Name),] # 16963

# Read in methylation file # by chromosome 
for(i in 1:22){ 
  meth=readRDS(paste0("GS20k_chr", i, "_mvals.rds"))
  
  # Subset to those with phenotype in question 
  meth=meth[,which(colnames(meth) %in% phenos$Sample_Sentrix_ID)] #80545 17464
  
  # Subset to probes passing QC 
  probes=read.table("/Cluster_Filespace/Marioni_Group/Elena/gs_osca/data/cpgs_tokeep.txt", header=F)
  meth=meth[which(row.names(meth) %in% probes$V1),] # 73108 17464
  
  # # Subset to 450k only 
  # meth=meth[which(row.names(meth)%in%anno$Name),]
  
  # Match order of IDs in phenotype and methylation file 
  ids=phenos$Sample_Sentrix_ID
  meth=meth[,match(ids,colnames(meth))]
  
  # Check order of IDs match between phenotype and methylation files 
  print(table(colnames(meth)==phenos$Sample_Sentrix_ID))
  
  # Convert to beta values 
  meth=m2beta(meth)
  
  # Mean impute - cannot have NAs in final file 
  meth <- apply(meth,1,meanimpute)
  
  # Transpose back original format - as apply changes the format
  meth=t(meth)
  
  # Prepare covariate matrix for regressions 
  # Match order of IDs with other files 
  samps <- samps_ref
  samps=samps[match(ids,samps$Sample_Sentrix_ID),] 
  
  # Check order of IDs match with other files 
  print(table(samps$Sample_Sentrix_ID==phenos$Sample_Sentrix_ID))
  print(table(samps$Sample_Sentrix_ID==colnames(meth)))
  
  ## Regression step - residualise for age, sex and batch 
  design.resid <- model.matrix(~sex + age + Batch, data=samps)
  fit.resid <- limma::lmFit(meth, design.resid)
  gc()
  meth <- limma::residuals.MArrayLM(fit.resid, meth)
  meth <- meth[!is.infinite(rowSums(meth)), ]
  rm(fit.resid)
  gc()
  
  ## Scale 
  meth=t(apply(meth,1,scale))
  
  ## Write out CpGs 
  cpgs=as.data.frame(row.names(meth))
  names(cpgs)[1]="CpG"
  fwrite(cpgs, paste0("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/cpgs_chr_index_bnp/GS20k_chr", i, "_subset_cpgs.txt"),row.names=F)
  
  # Save out residualised file 
  fwrite(meth, paste0("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/meth_resid_bnp/GS20k_chr", i, "_resid.txt"),row.names=F)
  
  ## Remove methylation object and clean up environment 
  rm(meth)
  gc()
  print(i)
} 
# BNP

## Write out ID order for later
saveRDS(ids, "/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/bnp_meth_ids.rds")


###############################################################
######## STEP 4 - PREPARE FIXED EFFECT COVARIATE FILES ########
###############################################################

#### HERE YOU MAY WANT TO TRUNCATE TO JUST WBCS - I JUST ADAPTED THIS TO MOST COVARIATES AS I HAD VARIOUS MODELS TO CHECK 

# Read in WBCs
wbc=read.table("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/wbc_quant.qcov", header=T)
names(wbc)[2] <- 'Sample_Sentrix_ID'

## GDF
ids <- readRDS("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/gdf_meth_ids.rds")
wbc1 <- wbc
gdf <- read.csv("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/GDF15_protein_file.csv")
wbc1 <- left_join(gdf, wbc1, by = 'Sample_Sentrix_ID')
wbc1=wbc1[match(ids,wbc1$FID.y),]
identical(wbc1$FID.y, gdf$Sample_Sentrix_ID)

## Remove Sample ID information 
# wbc1$FID=NULL
# wbc1$IID=NULL 
# wbc1$Sample_Name=NULL 
# wbc1$sex=ifelse(wbc1$sex%in%"M",0,1)
# wbc1$Batch=as.numeric(as.factor(wbc1$Batch))
wbc1 <- wbc1[c(35:39)]
wbc1[1:ncol(wbc1)]=apply(wbc1[,1:ncol(wbc1)],2,scale)

# Write out covariate file
write.table(x = as.matrix(wbc1),file = "/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/gdf_wbc_covariates.csv" ,quote = F, sep = ",", row.names = F, col.names = F)



# Read in WBCs
wbc=read.table("/Cluster_Filespace/Marioni_Group/Rob/EWAS_Disease_GS/wbc_quant.qcov", header=T)
names(wbc)[2] <- 'Sample_Sentrix_ID'

## BNP
ids <- readRDS("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/bnp_meth_ids.rds")
wbc1 <- wbc
bnp <- read.csv("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/BNP_protein_file.csv")
wbc1 <- left_join(bnp, wbc1, by = 'Sample_Sentrix_ID')
wbc1=wbc1[match(ids,wbc1$FID.y),]
identical(wbc1$FID.y, bnp$Sample_Sentrix_ID)

## Remove Sample ID information 
# wbc1$FID=NULL
# wbc1$IID=NULL 
# wbc1$Sample_Name=NULL 
# wbc1$sex=ifelse(wbc1$sex%in%"M",0,1)
# wbc1$Batch=as.numeric(as.factor(wbc1$Batch))
wbc1 <- wbc1[c(35:39)]
wbc1[1:ncol(wbc1)]=apply(wbc1[,1:ncol(wbc1)],2,scale)

# Write out covariate file
write.table(x = as.matrix(wbc1),file = "/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/bnp_wbc_covariates.csv" ,quote = F, sep = ",", row.names = F, col.names = F)

##############################################################
######## STEP 5 - COMBINE INDIVIDUAL CHROMOSOME FILES ########
##############################################################

## GDF

library(data.table)
library(limma)
library(lumi)
library(readxl)
library(tidyverse)

# Change working directory 
# setwd("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/meth_resid_gdf/") 
setwd('/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/meth_resid_gdf/')

# Extract files 
files=list.files(".", ".txt")
files=files[order(files)]

# Read in and rbind all methylation files 
data <- rbindlist(lapply(files,fread)) # 752722 17489

## Write out final file
# fwrite(x =as.matrix(data), "/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/Methylation/gdf_pred_17489.csv", sep = ",", row.names = F, col.names = F, quote = F) 
fwrite(x =as.matrix(data), "/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/Methylation_GDF/gdf_pred_17489.csv", sep = ",", row.names = F, col.names = F, quote = F) 


# Create list of cpgs 
setwd('/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/cpgs_chr_index_gdf/')

# Extract CpGs - will need this for processing final results files as the row.names get lost in methylation file in BayesR
cpgs=list.files("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/cpgs_chr_index_gdf/", ".txt")
cpgs=cpgs[order(cpgs)]

# Ensure that order is same as methylation file just created 
# methylation
ids.cpg=gsub("GS20k_chr", "", files)
ids.cpg=gsub("_.*", "", ids.cpg)
# cpgs
ids1.cpg=gsub("GS20k_chr", "", cpgs)
ids1.cpg=gsub("_.*", "", ids1.cpg)
# is order same? 
table(ids.cpg==ids1.cpg)

# Read in cpg lists 
# Change working directory 
# setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/CpGs/") 
cg <- rbindlist(lapply(cpgs,fread))
names(cg)[1]="Marker"


# Save out file
fwrite(x =as.matrix(cg), "/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/CpGs_GDF/gdf_cpg_subset_list.csv", sep = ",", row.names = F, col.names = T, quote = F) 




## BNP

library(data.table)
library(data.table)
library(limma)
library(lumi)
library(readxl)
library(tidyverse)

# Change working directory 
setwd('/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/meth_resid_bnp/')

# Extract files 
files=list.files(".", ".txt")
files=files[order(files)]

# Read in and rbind all methylation files 
data <- rbindlist(lapply(files,fread)) # 752722 16963

## Write out final file
fwrite(x =as.matrix(data), "/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/Methylation_BNP/bnp_pred_16.csv", sep = ",", row.names = F, col.names = F, quote = F) 



# Create list of cpgs 
setwd('/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/cpgs_chr_index_bnp/')

# Extract CpGs - will need this for processing final results files as the row.names get lost in methylation file in BayesR
cpgs=list.files("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/cpgs_chr_index_bnp/", ".txt")
cpgs=cpgs[order(cpgs)]

# Ensure that order is same as methylation file just created 
# methylation
ids.cpg=gsub("GS20k_chr", "", files)
ids.cpg=gsub("_.*", "", ids.cpg)
# cpgs
ids1.cpg=gsub("GS20k_chr", "", cpgs)
ids1.cpg=gsub("_.*", "", ids1.cpg)
# is order same? 
table(ids.cpg==ids1.cpg)

# Read in cpg lists 
# Change working directory 
# setwd("/Cluster_Filespace/Marioni_Group/Rob/CRP/BayesR/CpGs/") 
cg <- rbindlist(lapply(cpgs,fread))
names(cg)[1]="Marker"

# Save out file
fwrite(x =as.matrix(cg), "/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/CpGs_BNP/bnp_cpg_subset_list.csv", sep = ",", row.names = F, col.names = T, quote = F) 


#####################################
#### STEP 7 - BAYESR  ###############
#####################################

screen
 
## GDF
cd /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/
  
../../../BayesRRcmd/src/brr --data-file Methylation/gdf_pred_17489.csv --pheno Phenotypes/gdf_resid.csvphen --analysis-type preprocess --fixed_effects gdf_wbc_covariates.csv --fixedEffectNumber 5 --thread 12 --thread-spawned 12 --marker-cache --seed 1 

../../../BayesRRcmd/src/brr --data-file Methylation/gdf_pred_17489.csv --pheno Phenotypes/gdf_resid.csvphen --fixed_effects gdf_wbc_covariates.csv --fixedEffectNumber 5 --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --S "0.001,0.01,0.1" --mcmc-samples Outputs_gdf/gdf_resid.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1



screen

## BNP
cd /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/
  
../../../BayesRRcmd/src/brr --data-file Methylation/bnp_pred_16963.csv --pheno Phenotypes/bnp_resid.csvphen --analysis-type preprocess --fixed_effects bnp_wbc_covariates.csv --fixedEffectNumber 5 --thread 12 --thread-spawned 12 --marker-cache --seed 1 

../../../BayesRRcmd/src/brr --data-file Methylation/bnp_pred_16963.csv--pheno Phenotypes/bnp_resid.csvphen --fixed_effects bnp_wbc_covariates.csv --fixedEffectNumber 5 --analysis-type ppbayes --chain-length 10000 --burn-in 5000 --thin 5 --S "0.001,0.01,0.1" --mcmc-samples Outputs_bnp/bnp_resid.csv --thread 12 --thread-spawned 12 --marker-cache --seed 1

