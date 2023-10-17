####################################################################################

### Project EpiScores into W4 test set (GS)

####################################################################################

# Project episcores for both proteins into GS and save out for cox models
# Calculate performance metrics in GS (measured protein correlations, incremental R2)

####################################################################################

### W4 test set - EPIC array first, 10cv, pQTL unadjusted 

####################################################################################

cd /Local_Scratch/Danni/GDFBNP/

screen

R

library(tidyverse)

# Weights 10cv 
weights_gdf_10cv_pQTL_450k <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/gdf15W1W3_train_weights_450K_10cv_pQTL_unadjusted.csv")
weights_bnp_10cv_pQTL_450k <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/nt.probnpW1W3_train_weights_450K_10cv_pQTL_unadjusted.csv")

# Test sets EPIC
BNP_x_test <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/BNP_test_x_subset_EPIC.rds")
GDF_x_test <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/GDF_test_x_subset_EPIC.rds")

# Test sets 450K
BNP_x_test_450 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/BNP_test_x_subset.rds")
GDF_x_test_450 <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/GDF_test_x_subset.rds")

# Read in proteins that have not been regressed onto any covariates 
BNP_test <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W4_BNP_pQTL.csv")
GDF_test <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/W4_GDF15_pQTL.csv")

# Read in unrelated W4 individuals to training 
rel1 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/unrelated_W4_GDF.csv")
rel2 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/unrelated_W4_BNP.csv")

# Check no related individuals in test set 
which(!GDF_test$Sample_Name %in% rel1$Sample_Name) # 0
which(!BNP_test$Sample_Name %in% rel2$Sample_Name) # 0

# Read in PGS 
GDFPGS <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gs_genomewide_gdf.all.score")
BNPPGS <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gs_genomewide_bnp.all.score")

# Set location
location <- "/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores/"

##############################################################################################

### GDF15 - EPIC, 10cv

##############################################################################################

# Set y test variable
ytest <- GDF_test

# Set x test varuiabe
xtest <- GDF_x_test
xtest <- xtest[which(rownames(xtest) %in% ytest$Sample_Sentrix_ID),]

# Assign coefs
coefs <- weights_gdf_10cv_pQTL

  # GENERATE SCORES
  i <- 2
  names(ytest)[2] <- "pheno" # Assign a generic name to the protein variable

  xtest <- t(xtest) # transpose so rows are DNAm and people are columns 
  overlap <- which(rownames(xtest) %in% coefs$CpG) # Find the overlap between CpGs in the predictor weights column and the CpGs in the test methylation data 
  xtest <- xtest[overlap,] # Subset methylation CpG sites based on this overlap 
  match <- xtest[match(coefs$CpG, rownames(xtest)),] # Match up the order of CpGs in methylation file to those in the CpG predictor weights column
  calc <- match * coefs[,3] # Multiply beta predictor weights to the CpG methylation values for each person (column) in the methylation dataset 
  calc <- na.omit(calc)
  sum <- colSums(calc) # Sum the score for each person (column) in the methylation dataset 
  export_sum <- as.data.frame(sum) # Get the scores ready to be written to file
  names(export_sum)[1] <- "Scores"
  export_sum$Predictor <- "GDF15"
  export_sum$Sample_Name <- rownames(export_sum)


  # GENERATE INCREMENTAL R2
  target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
  names(export_sum)[3] <- "Sample_Sentrix_ID"
  join <- left_join(ytest, export_sum, by = "Sample_Sentrix_ID")
  join <- left_join(join, target, by = "Sample_Sentrix_ID")

  null <- summary(lm(pheno ~ age + sex, data=join))$r.squared
  full <- summary(lm(pheno ~ age + sex + Scores, data=join))$r.squared
  print(round(100*(full - null), 3))

# [1] 12.232

# Read in export sum file and join to protein y 
  names(export_sum)[3] <- "Sample_Sentrix_ID"
  join <- left_join(export_sum, ytest, by = "Sample_Sentrix_ID")
  join <- left_join(join, target, by = "Sample_Sentrix_ID")

# Now calculate correlations 
cor.test(join$pheno, join$Scores, method = "pearson")


# Plot correlation between protein and EpiScore
library(ggpubr)
plot1 <- ggplot(join, aes(x=pheno, y=Scores)) +
geom_point(colour = "red", size = 0.5) +
geom_smooth(method='lm', colour = "red") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) + ylim(-1,2) +
xlab("GDF15 protein") + ylab("GDF15 EpiScore EPIC") # , label.x = 4.5, label.y = 4.4

write.csv(export_sum, file = paste0(location, "projections_GS_GDF_scores_W4_EPIC_10cv.csv"), row.names = F)

# Join in PRS and try those in the models
GDFPGS <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gs_genomewide_gdf.all.score", header = T)
BNPPGS <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gs_genomewide_bnp.all.score", header = T)
names(GDFPGS)[2] <- 'Sample_Name'
names(BNPPGS)[2] <- 'Sample_Name'
GDFPGS <- GDFPGS[-1]
BNPPGS <- BNPPGS[-1]
names(GDFPGS)[2] <- 'GDFPRS'
names(BNPPGS)[2] <- 'BNPPRS'

join <- left_join(join, GDFPGS, by = 'Sample_Name')
join <- left_join(join, BNPPGS, by = 'Sample_Name')

null <- summary(lm(pheno ~ age + sex, data=join))$r.squared
null_PRS <- summary(lm(pheno ~ age + sex + GDFPRS, data=join))$r.squared
null_score <- summary(lm(pheno ~ age + sex + Scores, data=join))$r.squared
full <- summary(lm(pheno ~ age + sex + Scores + GDFPRS, data=join))$r.squared
  print(round(100*(full - null), 3)) # 15.77



##############################################################################################

### BNP - EPIC, 10cv

##############################################################################################

# Set y test variable
ytest <- BNP_test

# Set x test variable
xtest <- BNP_x_test
xtest <- xtest[which(rownames(xtest) %in% ytest$Sample_Sentrix_ID),]

# Assign coefs
coefs <- weights_bnp_10cv_pQTL

  # GENERATE SCORES
  i <- 1
  names(ytest)[2] <- "pheno" # Assign a generic name to the protein variable

  xtest <- t(xtest) # transpose so rows are DNAm and people are columns 
  overlap <- which(rownames(xtest) %in% coefs$CpG) # Find the overlap between CpGs in the predictor weights column and the CpGs in the test methylation data 
  xtest <- xtest[overlap,] # Subset methylation CpG sites based on this overlap 
  match <- xtest[match(coefs$CpG, rownames(xtest)),] # Match up the order of CpGs in methylation file to those in the CpG predictor weights column
  calc <- match * coefs[,3] # Multiply beta predictor weights to the CpG methylation values for each person (column) in the methylation dataset 
  calc <- na.omit(calc)
  sum <- colSums(calc) # Sum the score for each person (column) in the methylation dataset 
  export_sum <- as.data.frame(sum) # Get the scores ready to be written to file
  names(export_sum)[1] <- "Scores"
  export_sum$Predictor <- "GDF15"
  export_sum$Sample_Name <- rownames(export_sum)


  # GENERATE INCREMENTAL R2
  target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
  names(export_sum)[3] <- "Sample_Sentrix_ID"
  join <- left_join(ytest, export_sum, by = "Sample_Sentrix_ID")
  join <- left_join(join, target, by = "Sample_Sentrix_ID")
  null <- summary(lm(pheno ~ age + sex, data=join))$r.squared
  full <- summary(lm(pheno ~ age + sex + Scores, data=join))$r.squared
  print(round(100*(full - null), 3))

# 5.735

# Read in export sum file and join to protein y 
  names(export_sum)[3] <- "Sample_Sentrix_ID"
  join <- left_join(export_sum, ytest, by = "Sample_Sentrix_ID")
  join <- left_join(join, target, by = "Sample_Sentrix_ID")

# Now calculate correlations 
cor.test(join$pheno, join$Scores, method = "pearson")

# # Plot correlation between protein and EpiScore

plot3 <- ggplot(join, aes(x=pheno, y=Scores)) +
geom_point(colour = "red", size = 0.5) +
geom_smooth(method='lm', colour = "red") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) +
xlab("NT-proBNP protein") + ylab("NT-proBNP EpiScore EPIC") # , label.x = 4.5, label.y = 4.4


write.csv(export_sum, file = paste0(location, "projections_GS_BNP_scores_W4_EPIC_10cv.csv"), row.names = F)

# Join in PRS and try those in the models
GDFPGS <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gs_genomewide_gdf.all.score", header = T)
BNPPGS <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gs_genomewide_bnp.all.score", header = T)
names(GDFPGS)[2] <- 'Sample_Name'
names(BNPPGS)[2] <- 'Sample_Name'
GDFPGS <- GDFPGS[-1]
BNPPGS <- BNPPGS[-1]
names(GDFPGS)[2] <- 'GDFPRS'
names(BNPPGS)[2] <- 'BNPPRS'

join <- left_join(join, GDFPGS, by = 'Sample_Name')
join <- left_join(join, BNPPGS, by = 'Sample_Name')

null <- summary(lm(pheno ~ age + sex, data=join))$r.squared
null_PRS <- summary(lm(pheno ~ age + sex + BNPPRS, data=join))$r.squared
null_score <- summary(lm(pheno ~ age + sex + Scores, data=join))$r.squared
full <- summary(lm(pheno ~ age + sex + Scores + BNPPRS, data=join))$r.squared
  print(round(100*(full - null), 3)) # 6.94


##############################################################################################

### GDF15 - 450k, 10cv

##############################################################################################

# Set y test variable
ytest <- GDF_test

# Set x test varuiabe
xtest <- GDF_x_test_450
xtest <- xtest[which(rownames(xtest) %in% ytest$Sample_Sentrix_ID),]

# Assign coefs
coefs <- weights_gdf_10cv_pQTL_450k

  # GENERATE SCORES
  i <- 2
  names(ytest)[2] <- "pheno" # Assign a generic name to the protein variable

  xtest <- t(xtest) # transpose so rows are DNAm and people are columns 
  overlap <- which(rownames(xtest) %in% coefs$CpG) # Find the overlap between CpGs in the predictor weights column and the CpGs in the test methylation data 
  xtest <- xtest[overlap,] # Subset methylation CpG sites based on this overlap 
  match <- xtest[match(coefs$CpG, rownames(xtest)),] # Match up the order of CpGs in methylation file to those in the CpG predictor weights column
  calc <- match * coefs[,3] # Multiply beta predictor weights to the CpG methylation values for each person (column) in the methylation dataset 
  calc <- na.omit(calc)
  sum <- colSums(calc) # Sum the score for each person (column) in the methylation dataset 
  export_sum <- as.data.frame(sum) # Get the scores ready to be written to file
  names(export_sum)[1] <- "Scores"
  export_sum$Predictor <- "GDF15"
  export_sum$Sample_Name <- rownames(export_sum)


  # GENERATE INCREMENTAL R2
  target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
  names(export_sum)[3] <- "Sample_Sentrix_ID"
  join <- left_join(ytest, export_sum, by = "Sample_Sentrix_ID")
  join <- left_join(join, target, by = "Sample_Sentrix_ID")

  null <- summary(lm(pheno ~ age + sex, data=join))$r.squared
  full <- summary(lm(pheno ~ age + sex + Scores, data=join))$r.squared
  print(round(100*(full - null), 3))

# 1] 11.333

# Read in export sum file and join to protein y 
  names(export_sum)[3] <- "Sample_Sentrix_ID"
  join <- left_join(export_sum, ytest, by = "Sample_Sentrix_ID")
  join <- left_join(join, target, by = "Sample_Sentrix_ID")

# Now calculate correlations 
cor.test(join$pheno, join$Scores, method = "pearson")



# Plot correlation between protein and EpiScore
library(ggpubr)
plot2 <- ggplot(join, aes(x=pheno, y=Scores)) +
geom_point(colour = "tan1", size = 0.5) +
geom_smooth(method='lm', colour = "tan1") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) + ylim(-1,2) +
xlab("GDF15 protein") + ylab("GDF15 EpiScore 450k") # , label.x = 4.5, label.y = 4.4


write.csv(export_sum, file = paste0(location, "projections_GS_GDF_scores_W4_450k_10cv.csv"), row.names = F)

# Join in PRS and try those in the models
GDFPGS <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gs_genomewide_gdf.all.score", header = T)
BNPPGS <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gs_genomewide_bnp.all.score", header = T)
names(GDFPGS)[2] <- 'Sample_Name'
names(BNPPGS)[2] <- 'Sample_Name'
GDFPGS <- GDFPGS[-1]
BNPPGS <- BNPPGS[-1]
names(GDFPGS)[2] <- 'GDFPRS'
names(BNPPGS)[2] <- 'BNPPRS'

join <- left_join(join, GDFPGS, by = 'Sample_Name')
join <- left_join(join, BNPPGS, by = 'Sample_Name')

null <- summary(lm(pheno ~ age + sex, data=join))$r.squared
null_PRS <- summary(lm(pheno ~ age + sex + GDFPRS, data=join))$r.squared
null_score <- summary(lm(pheno ~ age + sex + Scores, data=join))$r.squared
full <- summary(lm(pheno ~ age + sex + Scores + GDFPRS, data=join))$r.squared
  print(round(100*(full - null), 3)) # 15.17


##############################################################################################

### BNP - 450k, 10cv

##############################################################################################

# Set y test variable
ytest <- BNP_test

# Set x test variable
xtest <- BNP_x_test_450
xtest <- xtest[which(rownames(xtest) %in% ytest$Sample_Sentrix_ID),]

# Assign coefs
coefs <- weights_bnp_10cv_pQTL_450k

  # GENERATE SCORES
  i <- 1
  names(ytest)[2] <- "pheno" # Assign a generic name to the protein variable

  xtest <- t(xtest) # transpose so rows are DNAm and people are columns 
  overlap <- which(rownames(xtest) %in% coefs$CpG) # Find the overlap between CpGs in the predictor weights column and the CpGs in the test methylation data 
  xtest <- xtest[overlap,] # Subset methylation CpG sites based on this overlap 
  match <- xtest[match(coefs$CpG, rownames(xtest)),] # Match up the order of CpGs in methylation file to those in the CpG predictor weights column
  calc <- match * coefs[,3] # Multiply beta predictor weights to the CpG methylation values for each person (column) in the methylation dataset 
  calc <- na.omit(calc)
  sum <- colSums(calc) # Sum the score for each person (column) in the methylation dataset 
  export_sum <- as.data.frame(sum) # Get the scores ready to be written to file
  names(export_sum)[1] <- "Scores"
  export_sum$Predictor <- "GDF15"
  export_sum$Sample_Name <- rownames(export_sum)


  # GENERATE INCREMENTAL R2
  target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
  names(export_sum)[3] <- "Sample_Sentrix_ID"
  join <- left_join(ytest, export_sum, by = "Sample_Sentrix_ID")
  join <- left_join(join, target, by = "Sample_Sentrix_ID")
  null <- summary(lm(pheno ~ age + sex, data=join))$r.squared
  full <- summary(lm(pheno ~ age + sex + Scores, data=join))$r.squared
  print(round(100*(full - null), 3))

# 5.37

# Read in export sum file and join to protein y 
  names(export_sum)[3] <- "Sample_Sentrix_ID"
  join <- left_join(export_sum, ytest, by = "Sample_Sentrix_ID")
  join <- left_join(join, target, by = "Sample_Sentrix_ID")

# Now calculate correlations 
cor.test(join$pheno, join$Scores, method = "pearson")



# # Plot correlation between protein and EpiScore

plot4 <- ggplot(join, aes(x=pheno, y=Scores)) +
geom_point(colour = "tan1", size = 0.5) +
geom_smooth(method='lm', colour = "tan1") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) +
xlab("NT-proBNP protein") + ylab("NT-proBNP EpiScore 450k") # , label.x = 4.5, label.y = 4.4


write.csv(export_sum, file = paste0(location, "projections_GS_BNP_scores_W4_450k_10cv.csv"), row.names = F)

# Join in PRS and try those in the models
GDFPGS <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gs_genomewide_gdf.all.score", header = T)
BNPPGS <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gs_genomewide_bnp.all.score", header = T)
names(GDFPGS)[2] <- 'Sample_Name'
names(BNPPGS)[2] <- 'Sample_Name'
GDFPGS <- GDFPGS[-1]
BNPPGS <- BNPPGS[-1]
names(GDFPGS)[2] <- 'GDFPRS'
names(BNPPGS)[2] <- 'BNPPRS'

join <- left_join(join, GDFPGS, by = 'Sample_Name')
join <- left_join(join, BNPPGS, by = 'Sample_Name')

null <- summary(lm(pheno ~ age + sex, data=join))$r.squared
null_PRS <- summary(lm(pheno ~ age + sex + BNPPRS, data=join))$r.squared
null_score <- summary(lm(pheno ~ age + sex + Scores, data=join))$r.squared
full <- summary(lm(pheno ~ age + sex + Scores + BNPPRS, data=join))$r.squared
  print(round(100*(full - null), 3)) # 7.093




#######################################################################

### Patchwork plots together for test correlations

library(patchwork)
location <- '/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores/'

pdf(paste0(location, "correlation_plots_GS_EPIC.pdf"), width = 14, height = 10)
plot1 + plot3 + plot2 + plot4
dev.off()


#######################################################################

### Process weights files for suppl tables

# Suppl table for within GS

gdf450 <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/gdf15W1W3_train_weights_450K_10cv_pQTL_unadjusted.csv")
gdfepic <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/gdf15W1W3_train_weights_EPIC_10cv_pQTL_unadjusted.csv")

bnp450 <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/nt.probnpW1W3_train_weights_450K_10cv_pQTL_unadjusted.csv")
bnpepic <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/nt.probnpW1W3_train_weights_EPIC_10cv_pQTL_unadjusted.csv")

dim(gdf450)
dim(bnp450)
dim(gdfepic)
dim(bnpepic)

gdf450$Array <- '450k'
bnp450$Array <- '450k'
gdfepic$Array <- 'EPIC'
bnpepic$Array <- 'EPIC'

bind <- rbind(gdf450, bnp450)
bind <- rbind(bind, gdfepic)
bind <- rbind(bind, bnpepic)

write.csv(bind, '/Local_Scratch/Danni/GDFBNP/00_Updated_results/suppl_table_within_GS_predictor_weights.csv', row.names = F)


# Suppl table for whole of GS

gdf450 <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/gdf1520k_450K_20cv_pQTL_unadjusted.csv")
gdfepic <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/gdf1520k_EPIC_20cv.csv")

bnp450 <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/nt.probnp20k_450K_20cv_pQTL_unadjusted.csv")
bnpepic <- read.csv("/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/nt.probnpW1W3_train_weights_EPIC.csv")

dim(gdf450)
dim(bnp450)
dim(gdfepic)
dim(bnpepic)


gdf450$Array <- '450k'
bnp450$Array <- '450k'
gdfepic$Array <- 'EPIC'
bnpepic$Array <- 'EPIC'

bind <- rbind(gdf450, bnp450)
bind <- rbind(bind, gdfepic)
bind <- rbind(bind, bnpepic)

write.csv(bind, '/Local_Scratch/Danni/GDFBNP/00_Updated_results/suppl_table_whole_GS_predictor_weights.csv', row.names = F)



