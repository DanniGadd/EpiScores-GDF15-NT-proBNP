
####################################################################################

# Covid EpiScores analyses 

####################################################################################

library(tidyverse)
library(readxl)
library(lme4)

# Load covid data from archie
t <- read_excel("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/2021-09-03 C19 cases Feb.xlsx")
t <- as.data.frame(t)

# Take a look at variables available in this file (including those without DNAm data)
length(unique(t$id)) # 1713 unique individuals 
table(t$covid) # 554 with diagnosis of covid 
table(t$smr) # 31 in hospital 
table(t$icu) # 6 in ICU 

e <- t %>% filter(t$covid == "1") 
names(e)[1] <- "Sample_Name"

# > dim(e)
# [1] 554   8

# Now add in age at covid 
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/C19_test_dates_25Oct2021.csv")
names(g)[1] <- "Sample_Name"

e2 <- merge(e, g, by.x = 'Sample_Name', all.x = TRUE)

hasAntiBD <- !is.na(e2$S3_AntiBD_Year)
hasSwab <- !is.na(e2$S3_SwabDate_Year)
hasLinkedTest <- !is.na(e2$Linked_TestDate)

# Remove individuals who reported having covid in CL1 but not in CL2
Mismatch <- (e2$Had_COVID > 0) & (!e2$S2_Had_COVID > 0)

# Replace NA with 0 in covidLife1And2Mismatch
Mismatch[is.na(Mismatch)] <- 0
e2 <- e2[!Mismatch, ] # only one person has been removed - 553 now 


# Read in appt table to extract all baseline appointment dates (not just those in CovidLife)
apptTable <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/2021-07-30_appt.csv")

apptToTargetTableIndex <- match(e2$Sample_Name, apptTable$id)
apptDateString <- apptTable[apptToTargetTableIndex, 'appt']

# Extract Date from first half of the timestamp
apptDate <- lapply(apptDateString, function(x) {as.Date(strsplit(x, ' ')[[1]][[1]], '%Y-%m-%d')})

# Use date from the following sources if available: linked test, antibd, swab. Else use an approximate date of 01/01/2021
covidDate <- lapply(1:nrow(e2), function(rowName) {
  row <- e2[rowName, ]
  if (!is.na(row$Linked_TestDate)) {
    as.Date(row$Linked_TestDate, '%Y-%m-%d')
  } else if (!is.na(row$S3_AntiBD_Year)) {
    as.Date(paste(row$S3_AntiBD_Year, row$S3_AntiBD_Month, row$S3_AntiBD_Day, sep = '/'), '%Y/%m/%d')
  } else if (!is.na(row$S3_SwabDate_Year)) {
    as.Date(paste(row$S3_SwabDate_Year, row$S3_SwabDate_Month, row$S3_SwabDate_Day, sep = '/'), '%Y/%m/%d')
  } else {
    as.Date('2021-01-01', '%Y-%m-%d')
  }
})

# Difference between appointment (baseline) date and covid date
covidApptDiff <- sapply(1:length(apptDate), function(i) {as.numeric(covidDate[[i]] - apptDate[[i]]) / 365})

######################################################

### Load the protein data in 20k

# Add maximal covariate info from protein file
prot <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/file_input_030821/GS20K_GDF15_NT_proBNP.PHE")
names(prot)[2] <- "Sample_Name"
prot <- prot[-c(3:6)]

# Add in protein data as predictors - pre prepped in N used for analyses
GDF15 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/GDF15_data_cox.csv")
NT <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/BNP_data_cox.csv")

GDF15 <- GDF15[c(2,4)]
NT <- NT[c(2,6)]

d1 <- left_join(prot, GDF15, by = 'Sample_Name')
d1 <- left_join(d1, NT, by = 'Sample_Name')

table(is.na(d1$gdf15))
table(is.na(d1$nt.probnp))


# Join file into those with covid 
d2 <- left_join(e2, d1, by = "Sample_Name")
dim(d2)

######################################################

# Add covid appt difference onto baseline age
d2$covidAge <- d2$age + covidApptDiff

# Assign difference as extra column
d2$covidDiff <- covidApptDiff

# filter to just smr cases 
d2_smr <- d2[which(d2$smr %in% "1"),] # 29 cases

# calculate mean difference 
mean <- mean(d2_smr$covidDiff, na.rm = T) # 11.864
sd <- sd(d2_smr$covidDiff) # 1.354

######################################################

# GDF

# d2_GDF <- d2[complete.cases(d2$GDF_score),]
d2_GDF <- d2[complete.cases(d2$gdf15),]

# > table(d2_GDF$smr)

#   0   1
# 463  28

# > dim(d2_GDF)
# [1] 491  51

library(glm2)

d2_GDF$smr <- as.factor(d2_GDF$smr)
results <- data.frame(episcore = "X", outcome = "X", n = "X", Beta = "X", SE = "X", p = "X")
markers <- colnames(d2_GDF)[c(48)]
outcome <- "hospitalised covid"

vars <- 'smr'

list3=list()
# Loop through predictors and input disease of interest as var variable 
for(j in c("gdf15")){ 
for(var in vars[1]){ 
  name <- 'gdf15'
    list3[[j]][[var]] <- summary(glm(smr ~ scale(d2_GDF[,name]) + scale(covidAge) + factor(sex), data = d2_GDF, 
     family = binomial))$coefficients[2,]
  }
} 

# combine 
l4 <- do.call('c', list3)
# Tidy df 
l4=lapply(split(l4,sub('.*\\.', '', names(l4))),function(x) do.call(rbind, x))
l4=as.data.frame(do.call("rbind",l4))
l4$pred=gsub("\\..*", "", row.names(l4))
l4$trait=gsub(".*\\.", "", row.names(l4))
# Convert to OR and conf int 
l4$HCI=exp(l4[,1]+(1.96*l4[,2]))
l4$LCI=exp(l4[,1]-(1.96*l4[,2]))
l4$Odds=exp(l4[,1])
names(l4)=c("Beta","SE","t","P","Predictor","Trait","HCI","LCI","Odds")


res1 <- l4

# BNP

# d2_BNP <- d2[complete.cases(d2$BNP_score),]
d2_BNP <- d2[complete.cases(d2$nt.probnp),]

# > table(d2_BNP$smr)
#   0   1
# 449  26

# > dim(d2_BNP)
# [1] 475  51

library(glm2)

d2_BNP$smr <- as.factor(d2_BNP$smr)
results2 <- data.frame(episcore = "X", outcome = "X", n = "X", Beta = "X", SE = "X", p = "X")
markers2 <- colnames(d2_BNP)[c(49)]
outcome2 <- "hospitalised covid"

vars <- 'smr'

list3=list()
# Loop through predictors and input disease of interest as var variable 
for(j in c("nt.probnp")){ 
for(var in vars[1]){ 
  name <- 'nt.probnp'
    list3[[j]][[var]] <- summary(glm(smr ~ scale(d2_BNP[,name]) + scale(covidAge) + factor(sex), data = d2_BNP, 
     family = binomial))$coefficients[2,]
  }
} 

# combine 
l4 <- do.call('c', list3)
# Tidy df 
l4=lapply(split(l4,sub('.*\\.', '', names(l4))),function(x) do.call(rbind, x))
l4=as.data.frame(do.call("rbind",l4))
l4$pred=gsub("\\..*", "", row.names(l4))
l4$trait=gsub(".*\\.", "", row.names(l4))
# Convert to OR and conf int 
l4$HCI=exp(l4[,1]+(1.96*l4[,2]))
l4$LCI=exp(l4[,1]-(1.96*l4[,2]))
l4$Odds=exp(l4[,1])
names(l4)=c("Beta","SE","t","P","Predictor","Trait","HCI","LCI","Odds")


# Bind and order results by P 
result <- rbind(res1, l4)

# Write off results file 
write.csv(result, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/03_covid/20k_result_glm_hospitalisations.csv", row.names = F)

##################################################################################
