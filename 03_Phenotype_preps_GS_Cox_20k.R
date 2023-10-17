#############################################################################################################

### Prep GDF15 and Nt-pro-BNP levels - 20k (for cox models)

#############################################################################################################

cd /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/00_Preps/

screen

R

library(readxl)
library(tidyverse)

# Read in new variable data from Riccardo (taken from GS folder on datastore)
d1 <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/file_input_030821/GS20K_GDF15_NT_proBNP.PHE")
names(d1)[2] <- "Sample_Name"

# Get complete pheno info 
table(is.na(d1$gdf15_rnk)) # 18414
table(is.na(d1$nt.probnp_rnk)) # 17863
GDF15 <- d1[complete.cases(d1$gdf15),]
NT <- d1[complete.cases(d1$nt.probnp),] 

# Read in 20k DNAm targets and add info
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
length(complete.cases(target$Sample_Sentrix_ID)) # 18413
GDF15 <- left_join(GDF15, target, by = "Sample_Name")
NT <- left_join(NT, target, by = "Sample_Name")

GDF15 <- GDF15[complete.cases(GDF15$Sample_Sentrix_ID),] # 17489
NT <- NT[complete.cases(NT$Sample_Sentrix_ID),] # 16963

# Transform and scale proteins 
library(bestNormalize)
for(i in colnames(GDF15)[4]){ 
  GDF15[,i]<- orderNorm(GDF15[,i])$x.t
}

for(i in colnames(NT)[6]){ 
  NT[,i]<- orderNorm(NT[,i])$x.t
}

GDF15[,4] <- scale(GDF15[,4])
NT[,6] <- scale(NT[,6])

# Save out versions for cox
write.csv(GDF15, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/GDF15_data_cox.csv", row.names = F)
write.csv(NT, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/BNP_data_cox.csv", row.names = F)


