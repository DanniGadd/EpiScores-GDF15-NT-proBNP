
####################################################################################

# Covid EpiScores analyses - long covid variable 

####################################################################################

library(tidyverse)
library(readxl)

# Load covid data from daniel 
t <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/GS_CL3_Samples_LongCovid.csv")

# Take a look at variables available in this file (including those without DNAm data)
length(unique(t$id)) # 3180 unique individuals 
table(t$S3_Had_COVID) # 347 with diagnosis of covid 
table(t$S3_SymptomLength1st) # 338 

# Filter by removing the NAs for symptom indications
d <- t %>% filter(t$S3_SymptomLength1st != "NA") # 338 remaining 

# Construct the binary variable of < 4 and > 4 weeks 
d$binary <- ifelse(d$S3_SymptomLengthAll == 1 | d$S3_SymptomLengthAll == 2, 0, 1)

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
d2 <- left_join(d, d1, by = "Sample_Name")
dim(d2)

# > dim(d2)
# [1] 338  30

# Look at the number of coded 1 vs 0 variables in subsets used - in the 338 with covid 

# GDF 
d2_GDF <- d2[complete.cases(d2$gdf15),]
table(d2_GDF$S3_SymptomLengthAll)

#   1   2   3   4
# 111  75  47  42

table(d2_GDF$binary)

#   0   1
# 186  89

# BNP 
d2_BNP <- d2[complete.cases(d2$nt.probnp),]
table(d2_BNP$S3_SymptomLengthAll)

#   1   2   3   4
# 107  70  47  40

table(d2_BNP$binary)

#   0   1
# 177  87


#################################

# GDF 

### Add in the covidage calculation 
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/C19_test_dates_25Oct2021.csv")
names(g)[1] <- "Sample_Name"

e <- d2_GDF

e <- e[complete.cases(e$gdf15),] # 100 181

e2 <- merge(e, g, by = 'Sample_Name', all.x = TRUE)

hasAntiBD <- !is.na(e2$S3_AntiBD_Year)
hasSwab <- !is.na(e2$S3_SwabDate_Year)
hasLinkedTest <- !is.na(e2$Linked_TestDate)

# Remove individuals who reported having covid in CL1 but not in CL2
Mismatch <- (e2$Had_COVID > 0) & (!e2$S2_Had_COVID > 0)

# Replace NA with 0 in covidLife1And2Mismatch
Mismatch[is.na(Mismatch)] <- 0
e2 <- e2[!Mismatch, ] # three removed - now 117 and 56 logn covid 

table(e2$binary)

#   0   1
# 182  87


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

# Add covid appt difference onto baseline age
e2$covidAge <- e2$age + covidApptDiff

# Assign difference as extra column
e2$covidDiff <- covidApptDiff

# filter to just smr cases 
e2_smr <- e2[which(e2$binary %in% "1"),]

# calculate mean difference 
mean <- mean(e2_smr$covidDiff, na.rm = T) # 11.29
sd <- sd(e2_smr$covidDiff) # 1.24


# Run glm() models with the binary variable
library(glm2)


e2$binary <- as.factor(e2$binary)
results <- data.frame(prot = "X", outcome = "X", n = "X", Beta = "X", SE = "X", p = "X")
markers <- c("gdf15")
outcome <- "long covid"


vars <- 'binary'

list3=list()
# Loop through predictors and input disease of interest as var variable 
for(j in c("gdf15")){ 
for(var in vars[1]){ 
  name <- 'gdf15'
    list3[[j]][[var]] <- summary(glm(binary ~ scale(e2[,name]) + scale(covidAge) + factor(sex), data = e2, 
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

result_GDF <- l4


############################

# BNP

### Add in the covidage calculation 
g <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/LBC_proteins_Jan2021/COVID_update/C19_test_dates_25Oct2021.csv")
names(g)[1] <- "Sample_Name"

e <- d2_BNP

e <- e[complete.cases(e$BNP_score),]
e <- e[complete.cases(e$nt.probnp),] #

e2 <- merge(e, g, by = 'Sample_Name', all.x = TRUE)

hasAntiBD <- !is.na(e2$S3_AntiBD_Year)
hasSwab <- !is.na(e2$S3_SwabDate_Year)
hasLinkedTest <- !is.na(e2$Linked_TestDate)

# Remove individuals who reported having covid in CL1 but not in CL2
Mismatch <- (e2$Had_COVID > 0) & (!e2$S2_Had_COVID > 0)

# Replace NA with 0 in covidLife1And2Mismatch
Mismatch[is.na(Mismatch)] <- 0
e2 <- e2[!Mismatch, ] # three removed - now 117 and 56 logn covid 

table(e2$binary)

#   0   1
# 174  85


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

# Add covid appt difference onto baseline age
e2$covidAge <- e2$age + covidApptDiff

# Assign difference as extra column
e2$covidDiff <- covidApptDiff

# filter to just smr cases 
e2_smr <- e2[which(e2$binary %in% "1"),] 

# calculate mean difference 
mean <- mean(e2_smr$covidDiff, na.rm = T) 
sd <- sd(e2_smr$covidDiff) 


# Run glm() models with the binary variable
library(glm2)


e2$binary <- as.factor(e2$binary)
results <- data.frame(prot = "X", outcome = "X", n = "X", Beta = "X", SE = "X", p = "X")
markers <- c("nt.probnp")
outcome <- "long covid"


vars <- 'binary'

list3=list()
# Loop through predictors and input disease of interest as var variable 
for(j in c("nt.probnp")){ 
for(var in vars[1]){ 
  name <- 'nt.probnp'
    list3[[j]][[var]] <- summary(glm(binary ~ scale(e2[,name]) + scale(covidAge) + factor(sex), data = e2, 
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

result_BNP <- l4


# Bind and order results by P 
results <- rbind(result_GDF, result_BNP)
result <- results[order(results$P),]

# Write off results file 
write.csv(result, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/03_covid/20k_result_glm_long_covid.csv", row.names = F)


### COMBINE THE FILES 

hosp <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/03_covid/20k_result_glm_hospitalisations.csv")

long <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/03_covid/20k_result_glm_long_covid.csv")

library(tidyverse)
join <- rbind(hosp, long)
join <- join[order(join$P),]

# FDR 
join$FDR <- p.adjust(join$P, method = "BH")

write.csv(join, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/03_covid/results_joint_all_covid.csv")

