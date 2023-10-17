############################################################################################
############################################################################################
################# Projection into LBC1936 ##################################################
############################################################################################
############################################################################################

## Projections in LBC1936 W1 

library(tidyverse)

# Read in LBC data file
load("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/Beta_3525_norm_bgcorrect_0.001BetaThreshold_probefilter.RObject")

# Read in EPIC annotation file
anno <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/EPIC_AnnotationObject_df.rds")

# Subset annotation file to probes common to 450k and EPIC array
common_anno <- anno[which(anno$Methyl450_Loci == "TRUE"),]

# Subset methylation data file to those in the common annotation file 
dat1 <- dat[rownames(dat) %in% rownames(common_anno),] # cpgs are currently rows, people are cols 

# Read in target file with matched IDs for methylation data 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")

# Filter to wave 1
dat3 <- target %>% filter(WAVE == "1") 

# Filter to LBC36
dat3  <- dat3 %>% filter(cohort == "LBC36") 


# Subset DNAm data to these individuals from W1 of LBC36 only 
dat1 <- dat1[,which(colnames(dat1) %in% dat3$Basename)]


# Transpose data so CpGs are columns 
dat1 <- t(dat1)

# Replace NA values in methylation data with imputed means for tidyness (means for each CpG column imputed)
for(t in 1:ncol(dat1)){
     dat1[is.na(dat1[,t]), t] <- mean(dat1[,t], na.rm = TRUE)
}

# Assign variable for data to the LBC DNAm dataset 
data = dat1

# Transpose back so CpGs are rows again 
data <- t(data)

library(readxl)

## Start to Process Files 

message("1. Loading data") 

message("1.1 Loading Methylation data - rows to be CpGs and columns to be individuals") 

# Read in weights
GDF <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/Results_241122/gdf1520k_450K_20cv_pQTL_unadjusted.csv")
BNP <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/Results_241122/nt.probnp20k_450K_20cv_pQTL_unadjusted.csv")

cpgs <- rbind(GDF, BNP)
cpgs$Mean_Beta_Value <- 0
cpgs <- cpgs[c(2,3,4,1)]
names(cpgs)[1] <- "CpG_Site"
names(cpgs)[4] <- "Predictor"

## Check if Data needs to be Transposed

message("2. Quality Control and data Preparation") 

message("2.1 Checking if Row Names are CpG Sites") 

if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data<-t(data) 
}

message("2.2 Subsetting CpG sites to those required for Predictor Calculation") 

## Subset CpG sites to those present on list for predictors 

coef=data[intersect(rownames(data), cpgs$CpG_Site),]

## Check if Beta or M Values 

m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

coef<-if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef) } else { message("Suspect that Beta Values are present");coef}
# coef_age<-if((range(coef_age,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef_age) } else { message("Suspect that Beta Values are present");coef_age}

## Scale Data if Needed 

ids = colnames(coef)
scaled <- apply(coef, 1, function(x) sd(x,na.rm = T)) 

coef <-  if(range(scaled)[1] == 1 & range(scaled)[2] == 1) { 
    coef
  } else { 
    coef_scale <- apply(coef, 1, scale)
    coef_scale <- t(coef_scale)
    coef_scale <- as.data.frame(coef_scale)
    colnames(coef_scale) <- ids
   coef_scale
  } 


message("2.3 Find CpGs not present in uploaded file, add these with mean Beta Value for CpG site from Training Sample") 

## Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site

coef <- if(nrow(coef) == 2111) { message("All sites present"); coef } else if(nrow(coef)==0){ 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
} else { 
  missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)),c("CpG_Site","Mean_Beta_Value")]
  message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
  mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef))
  row.names(mat) <- unique(missing_cpgs$CpG_Site)
  colnames(mat) <- colnames(coef) 
  mat[is.na(mat)] <- 1
  missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
    missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
  } else {missing_cpgs
  }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat=mat*missing_cpgs1$Mean_Beta_Value
  coef=rbind(coef,mat) } 

message("2.4 Convert NA Values to Mean for each Probe") 

## Convert NAs to Mean Value for all individuals across each probe 

na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
  return(methyl)
}

coef <- t(apply(coef,1,function(x) na_to_mean(x)))
# coef_age <- t(apply(coef_age,1,function(x) na_to_mean(x)))


message("3. Calculating the Predictors") 

loop = unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = tmp2[,1]
  }
} 


out$ID <- row.names(out) 


## Save File and Finish Up 
message("Analysis Finished! Thank you for using our application. Output File is called \"out\"") 

## Save file 

out2 <- as.data.frame(out)

# write.csv(out, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Model_projections/LBC1921_162_projections.csv", row.names = F)
write.csv(out2, "/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores/projections_LBC_W1_241122.csv", row.names = F)


############################################################################################

## Projections in LBC1936 W4

# Subset methylation data file to those in the common annotation file 
dat1 <- dat[rownames(dat) %in% rownames(common_anno),] # cpgs are currently rows, people are cols 

# Read in target file with matched IDs for methylation data 
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/target_QC_age_sex_date.csv")

# Filter to wave 1
dat3 <- target %>% filter(WAVE == "4")

# Filter to LBC36
dat3  <- dat3 %>% filter(cohort == "LBC36")

# Subset DNAm data to these individuals from W1 of LBC36 only 
dat1 <- dat1[,which(colnames(dat1) %in% dat3$Basename)]

# Transpose data so CpGs are columns 
dat1 <- t(dat1)

# Replace NA values in methylation data with imputed means for tidyness (means for each CpG column imputed)
for(t in 1:ncol(dat1)){
     dat1[is.na(dat1[,t]), t] <- mean(dat1[,t], na.rm = TRUE)
}

# Assign variable for data to the LBC DNAm dataset 
data = dat1

# Transpose back so CpGs are rows again 
data <- t(data)

library(readxl)

## Start to Process Files 

message("1. Loading data") 

message("1.1 Loading Methylation data - rows to be CpGs and columns to be individuals") 

# Read in weights
GDF <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/gdf1520k_450K_20cv_pQTL_unadjusted.csv")
BNP <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/nt.probnp20k_450K_20cv_pQTL_unadjusted.csv")

cpgs <- rbind(GDF, BNP)
cpgs$Mean_Beta_Value <- 0
cpgs <- cpgs[c(2,3,4,1)]
names(cpgs)[1] <- "CpG_Site"
names(cpgs)[4] <- "Predictor"

## Check if Data needs to be Transposed

message("2. Quality Control and data Preparation") 

message("2.1 Checking if Row Names are CpG Sites") 

if(ncol(data) > nrow(data)){
  message("It seems that individuals are rows - data will be transposed!")
  data<-t(data) 
}

message("2.2 Subsetting CpG sites to those required for Predictor Calculation") 

## Subset CpG sites to those present on list for predictors 

coef=data[intersect(rownames(data), cpgs$CpG_Site),]

## Check if Beta or M Values 

m_to_beta <- function (val) 
{
  beta <- 2^val/(2^val + 1)
  return(beta)
}

coef<-if((range(coef,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef) } else { message("Suspect that Beta Values are present");coef}
# coef_age<-if((range(coef_age,na.rm=T)> 1)[[2]] == "TRUE") { message("Suspect that M Values are present. Converting to Beta Values");m_to_beta(coef_age) } else { message("Suspect that Beta Values are present");coef_age}

## Scale Data if Needed 

ids = colnames(coef)
scaled <- apply(coef, 1, function(x) sd(x,na.rm = T)) 

coef <-  if(range(scaled)[1] == 1 & range(scaled)[2] == 1) { 
    coef
  } else { 
    coef_scale <- apply(coef, 1, scale)
    coef_scale <- t(coef_scale)
    coef_scale <- as.data.frame(coef_scale)
    colnames(coef_scale) <- ids
   coef_scale
  } 


message("2.3 Find CpGs not present in uploaded file, add these with mean Beta Value for CpG site from Training Sample") 

## Identify CpGs missing from input dataframe, include them and provide values as mean methylation value at that site

coef <- if(nrow(coef) == 2111) { message("All sites present"); coef } else if(nrow(coef)==0){ 
  message("There Are No Necessary CpGs in The dataset - All Individuals Would Have Same Values For Predictors. Analysis Is Not Informative!")
} else { 
  missing_cpgs = cpgs[-which(cpgs$CpG_Site %in% rownames(coef)),c("CpG_Site","Mean_Beta_Value")]
  message(paste(length(unique(missing_cpgs$CpG_Site)), "unique sites are missing - add to dataset with mean Beta Value from Training Sample", sep = " "))
  mat = matrix(nrow=length(unique(missing_cpgs$CpG_Site)),ncol = ncol(coef))
  row.names(mat) <- unique(missing_cpgs$CpG_Site)
  colnames(mat) <- colnames(coef) 
  mat[is.na(mat)] <- 1
  missing_cpgs1 <- if(length(which(duplicated(missing_cpgs$CpG_Site))) > 1) { 
    missing_cpgs[-which(duplicated(missing_cpgs$CpG_Site)),]
  } else {missing_cpgs
  }  
  ids = unique(row.names(mat))
  missing_cpgs1 = missing_cpgs1[match(ids,missing_cpgs1$CpG_Site),]
  mat=mat*missing_cpgs1$Mean_Beta_Value
  coef=rbind(coef,mat) } 

message("2.4 Convert NA Values to Mean for each Probe") 

## Convert NAs to Mean Value for all individuals across each probe 

na_to_mean <-function(methyl){
  methyl[is.na(methyl)]<-mean(methyl,na.rm=T)
  return(methyl)
}

coef <- t(apply(coef,1,function(x) na_to_mean(x)))
# coef_age <- t(apply(coef_age,1,function(x) na_to_mean(x)))


message("3. Calculating the Predictors") 

loop = unique(cpgs$Predictor)
out <- data.frame()
for(i in loop){ 
  tmp=coef[intersect(row.names(coef),cpgs[cpgs$Predictor %in% i,"CpG_Site"]),]
  tmp_coef = cpgs[cpgs$Predictor %in% i, ]
  if(nrow(tmp_coef) > 1) { 
    tmp_coef = tmp_coef[match(row.names(tmp),tmp_coef$CpG_Site),]
    out[colnames(coef),i]=colSums(tmp_coef$Coefficient*tmp)
  } else {
    tmp2 = as.matrix(tmp)*tmp_coef$Coefficient 
    out[colnames(coef),i] = tmp2[,1]
  }
} 


out$ID <- row.names(out) 


## Save File and Finish Up 
message("Analysis Finished! Thank you for using our application. Output File is called \"out\"") 

## Save file 

out2 <- as.data.frame(out)

# write.csv(out, "/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/Model_projections/LBC1921_162_projections.csv", row.names = F)
write.csv(out2, "/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores/projections_LBC_W4_241122.csv", row.names = F)

