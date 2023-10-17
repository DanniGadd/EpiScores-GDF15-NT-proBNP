##################################################################################

### DEMENTIA INTEGRATION - INCIDENT LBC

##################################################################################

screen

R

# Load in packages 
library(survival)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)
library(haven)
library(foreign)
library(haven)


# Load dementia diagnoses from Paul 
data <- read.spss("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/LBC1936_DNAmPredictorsOfSerumGDF15_and_NTproBNP_DG_23NOV2022.sav", to.data.frame=TRUE)
dem <- read_sav("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/LBC1936_DNAmPredictorsOfSerumGDF15_and_NTproBNP_DG_DRAFT_DEMENTIA_ASCERTAINMENT_23NOV2022.sav")
data <- as.data.frame(data)
dem <- as.data.frame(dem)

# How many cases of dementia
table(dem$ConsensusDementiaDiagnosis)
table(dem$dementia_code)

#   0   1 
# 747 118

# Do we have dates for onset for all 118 
length(which(dem$age_at_dementia_diagnosis > 1)) # 120
sub <- dem[which(dem$age_at_dementia_diagnosis > 1),]
write.csv(sub, '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/dementia_LBC.csv', row.names = F) # 2 individuals with MCI here, rest dementia

# How many MCI/dem diagnoses do we have 
sub2 <- dem[which(dem$ConsensusDementiaDiagnosis > 1),] 
write.csv(sub2, '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/diagnoses_LBC.csv', row.names = F) # 2 individuals with MCI here, rest dementia

table(sub2$ConsensusDementiaDiagnosis)

###########################################

### FILTER TO W2 individuals as a baseline (W1 didnt consent to dementia followup)

# Filter to W2 LBC1936 target
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target[which(target$WAVE %in% '2'),]
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 
dat3 <- dat3 %>% select('Basename', 'set', 'ID')

# Subset dementia file to W2 individuals
sub <- dem[which(dem$lbc36no %in% dat3$ID),] # 801

# Look at dementia diagnoses in the context of W2
table(sub$ConsensusDementiaDiagnosis)
table(sub$dementia_code)

# Remove individuals with possible dementia or MCI 
sub <- sub[-which(sub$ConsensusDementiaDiagnosis %in% c('MCI', 'Possible dementia', '-999')),] # 780
table(sub$ConsensusDementiaDiagnosis)
table(sub$dementia_code)

names(dat3)[3] <- 'lbc36no'

library(haven)
new <- read_sav("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/original_data/DemAcertainmentScreeningAge.sav")

####################################################################################

### INCDIENT ANALYSES - W2 baseline

####################################################################################

#### PREP PHENOTYPE FILE WITH AGE ALIVE AND AGE DEATH INFO 

# Calculate age at cesor date for all individuals (dead or alive)
d1 <- sub
d1$YOB <- 1936
d1$MOB <- 07

d1 <- left_join(d1, new, by = 'lbc36no')
d1$age_censor <- d1$agedays_DemAsc_w6
d1$age_censor <- d1$age_censor / 365.25

# Map those that died before the censor date (i.e. age dead occurs before age at censor)
dead <- d1[complete.cases(d1$ageyrs_death),] # 290 died overall
dead$died_pre_censor <- ifelse(dead$ageyrs_death < dead$age_censor, 1, 0) 
table(dead$died_pre_censor) # all 290 died pre censor date 

# Assign dead status to main file
d1$dead <- ifelse(d1$ageyrs_death > 0, 1, 0)
d1$dead[is.na(d1$dead)] <- 0
table(d1$dead)
#   0   1
# 490 290

# Create variable with age dead if died, or age at censor if alive at censor 
d1$dead <- as.character(d1$dead)
d1$lbc36no <- as.character(d1$lbc36no)
d1$aged <- ifelse(d1$dead == '1', d1$ageyrs_death, d1$age_censor)

# Calculate age at baseline (W2) in years 
d1$age_baseline = d1$agedays_w2/365.25

# Get dementia diagnosis dates and join into a variable as the event 
dem <- sub[which(sub$dementia_code %in% '1'),] # 108
dem <- dem %>% select('lbc36no', 'age_at_dementia_diagnosis') # age in years when dementia first occured - complete for all 108 cases
names(dem)[2] <- 'age_diagnosis'
d1 <- left_join(d1, dem, by = 'lbc36no')

# Calculate time to event
cox <- d1
cox$age_at_event = ifelse(cox$dementia_code %in% 1, cox$age_diagnosis, (ifelse(cox$dead %in% 1 & cox$dementia_code %in% 0, cox$ageyrs_death, cox$aged)))
cox$tte = cox$age_at_event - cox$age_baseline
cox$tte = as.numeric(cox$tte)

# Check those with diagnoses in the uncertain zone 1 year pre/post baseline - there arent that many so will not affect things much 
length(which(cox$tte < -1))
# [1] 1
length(which(cox$tte < 0))
# [1] 2
length(which(cox$tte < 1))
# [1] 8

# Remove those with negative tte or tte between -1 and 0
# cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
cox$tte = ifelse(cox$tte < 0, 'NA', cox$tte)
cox$dementia_code = as.numeric(cox$dementia_code)
cox$tte<-as.numeric(cox$tte)

# Filter by age (not needed in LBC but good to remember)
cox = cox[cox$age_at_event >=65,] # 780 - all above 

# Add episcore markers that have been pre-corrected for DNAm technical variables and rank-inverse based transformed
scores <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/original_data/W2_LBC_episcores_adjusted_technical.csv")
names(scores)[3] <- 'Basename'

names(dat3)[3] <- 'lbc36no'
cox <- left_join(cox, dat3, by = 'lbc36no')
cox <- left_join(cox, scores, by = 'Basename')

# Check complete episcores
table(complete.cases(cox$GDF15_score)) # 775
table(complete.cases(cox$Nt.proBNP_score)) # 775

# Check GDF15 available at W2 in the dataset
table(complete.cases(cox$GDF15_W2)) # 744 true, 36 false 

# Check events
table(complete.cases(cox$dementia_code)) # 780


# Run basic models for each protein
mod1 = coxph(Surv(tte, dementia_code) ~ scale(GDF15_score) + factor(sex) + age_baseline, data = cox)
mod2 = coxph(Surv(tte, dementia_code) ~ scale(Nt.proBNP_score) + factor(sex) + age_baseline, data = cox)

summary(mod1)
summary(mod2)

# Mean tte for cases
cases <- cox[which(cox$dementia_code %in% 1),]
mean(cases$tte, na.rm = T)
sd(cases$tte, na.rm = T)

# Max tte across study Cox PH
max(cox$tte, na.rm = T)
