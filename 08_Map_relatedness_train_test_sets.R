####################################################################################

### RELATEDNESS ASSESSMENT GS

####################################################################################

# W4 case who has relative in W1/W3, we remove relatives from training 
# W4 control who has relative in W1/W3, then we remove the W4 control 

####################################################################################

# To identify W4 cases, we will first run cox models in W4 with the traits, to save case tables and index relatedness

screen

R

# Load in packages 
library(survival)
library(kinship2)
library(coxme)
library(readxl)
library(tidyverse)
library(gdata)

# Read in the prepped pedigree file to cluster and create kinship matrix
ped <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree_formatted.csv")
kin <- with(ped, pedigree(volid, father, mother, sex, famid=famid))
kin_model <- kinship(kin) 

# Load function to Extract Lmekin Results 
extract_coxme_table <- function (mod){
  beta <- mod$coefficients #$fixed is not needed
  nvar <- length(beta)
  nfrail <- nrow(mod$var) - nvar
  se <- sqrt(diag(mod$var)[nfrail + 1:nvar])
  z<- round(beta/se, 2)
  p<- 1 - pchisq((beta/se)^2, 1)
  table=data.frame(cbind(beta,se,z,p))
  return(table)
}


####################################################################################

#### PREP PHENOTYPE FILE WITH AGE ALIVE AND AGE DEATH INFO 

####################################################################################

## DEMOGRAPHICS

## Read in GS 20k age/sex base file and add prevalent data to it 
all <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/agemonths.csv")
names(all)[2] <- "age"
names(all)[1] <- "Sample_Name"
dim(all) # 24088   3
prevalent <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/disease.csv")
names(prevalent)[1] <- "Sample_Name"
merged_prev <- merge(all, prevalent, by = "Sample_Name", all = T)
dim(merged_prev) # 24092   115

## SURVIVAL INFO

# Read in deaths data and count how many individuals have died since oct 2020 (when GP data was last sampled)
age_dead <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Diabetes_update_codes_yipeng/2022-03-10_age_at_death.csv")
length(which(age_dead$dod_ym > 202010)) # 227 individuals since oct 2020 

## extract year/month of death as separate variables ##
age_dead$y_dead <- as.numeric(substr(age_dead$dod_ym, 1, 4))
age_dead$m_dead <- as.numeric(substr(age_dead$dod_ym, 5, 6))

# Assign the dead individuals as a separate subset
age_dead_include <- age_dead[which(age_dead$dod_ym <= 202010),] # 1350
age_dead_exclude <- age_dead[which(age_dead$dod_ym > 202010),] # 227

# Calculate a more exact estimate (by year and month) for age of death in the included 1350 individuals that died 
age_dead_include$y_diff <- age_dead_include$y_dead - age_dead_include$yob
age_dead_include$m_diff <- (age_dead_include$m_dead - age_dead_include$mob)/12
age_dead_include$diff <- age_dead_include$y_diff + age_dead_include$m_diff

# Work out those who are alive (i.e. not in the list of 1350 dead people from above)
age_alive <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_COVID/Diabetes_update_codes_yipeng/2022-03-10_age_alive.csv")
age_alive <- age_alive[c(1:3,6,8)]
age_alive = age_alive[!age_alive$id %in% age_dead_include$id,] # 22738 individuals who are not in the dead people we include - this covers the remainder of GS who did not die

# Ensure that all individuals in the 22738 sample are coded as alive and all individuals in the dead file are coded as such 
age_alive$dead <- 0
age_dead_include$dead <- 1

# Find age the 'alive' people were in oct 2020
age_alive$y_diff <- 2020 - age_alive$yob
age_alive$m_diff <- (10 - age_alive$mob)/12
age_alive$diff <- age_alive$y_diff + age_alive$m_diff

# The included dead people will have their age at death taken forward (i.e. pre 2020) as the 'aged' column in cox loops - this is now calculated
# The excluded dead people will be classed as alive and will have their age at 2020 taken forward as the 'aged' column - we have just calculated this as part of the wider group of alive individuals

# Subset to just the cols needed for joining dead and alive data
age_alive <- age_alive[c(1:4,8)] 
age_dead_include <- age_dead_include[c(1,2,3,7,13)]

# Bind the rows of the alive and dead files together for the whole GS sample 
names(age_dead_include) <- c("id", "yob", "mob", "dead", "aged")
names(age_alive) <- c("id", "yob", "mob", "dead", "aged")
age = rbind(age_alive, age_dead_include)
dim(age) # 24088     5
table(age$dead)

names(age)[1] <- "Sample_Name"

## Add survival info to the 20k base file 
d1 <- left_join(merged_prev, age, by = "Sample_Name")
dim(d1) # [1] 24092   119

# Create a subset with DOBMOB to use to filter cases by in the diabetes code processing below
d2 <- d1[c(1,2,116,117)]


####################################################################################

## Read in phenotypes to covary for in the fully-adjusted models 

# Load in the original d1 files

# Alcohol (units) and (usual)
alcohol <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/alcohol.csv")
names(alcohol)[1] <- "Sample_Name"

# BMI at baseline 
BMI <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/body.csv")
names(BMI)[1] <- "Sample_Name"

# SIMD
simd <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/clinical/SIMD.csv")
names(simd)[1] <- "Sample_Name"

# EA
ea <- read.csv("/Cluster_Filespace/Marioni_Group/GS/GS_dataset/PCQ/education.csv")
names(ea)[1] <- "Sample_Name"

# Join up the covariate data to the new d1 file in 20k 
d1 <- left_join(d1, alcohol, by = "Sample_Name")
d1 <- left_join(d1, BMI, by = "Sample_Name")
d1 <- left_join(d1, simd, by = "Sample_Name")
d1 <- left_join(d1, ea, by = "Sample_Name")

# Read in epismoekr 
w1 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave1_epismoker.rds")
w3 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave3_epismoker.rds")
w4 <- readRDS("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/DNAm_preps/wave4_epismoker.rds")
bind <- rbind(w1, w3)
bind <- rbind(bind, w4) # 18779 individuals with DNAm epismoker calculated 
bind$Sample_Sentrix_ID <- row.names(bind)
bind <- bind[-1]

# Join Sample_Name info 
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
bind <- left_join(bind, target, by = "Sample_Sentrix_ID")
bind <- bind[c(1,3)]

# Join to the data 
d1 <- left_join(d1, bind, by = "Sample_Name")

# ## Remove outlying values from key covairates used in models (>3.5 SDs from mean in each direction)
# list <- c("units", "usual", "bmi", "rank", "years")
# for(i in list){ 
  
#   cutoff1 = mean(d1[,i], na.rm = T) + 3.5*sd(d1[,i], na.rm = T)
#   cutoff2 = mean(d1[,i], na.rm = T) - 3.5*sd(d1[,i], na.rm = T)
  
#   d1[,i][which(d1[,i] > cutoff1 | d1[,i] < cutoff2)] <- NA 
# } 

####################################################################################

# Add maximal covariate info from protein file
prot <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/file_input_030821/GS20K_GDF15_NT_proBNP.PHE")
names(prot)[2] <- "Sample_Name"
prot <- prot[-c(3:6)]
d1 <- left_join(d1, prot, by = "Sample_Name")

# Add in protein data as predictors - pre prepped in N used for analyses
GDF15 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/GDF15_data_cox.csv")
NT <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/BNP_data_cox.csv")

GDF15 <- GDF15[c(2,4)]
NT <- NT[c(2,6)]

d1 <- left_join(d1, GDF15, by = 'Sample_Name')
d1 <- left_join(d1, NT, by = 'Sample_Name')

table(is.na(d1$gdf15))
table(is.na(d1$nt.probnp))

# Set the clock for variables of interest - the usual clock for running all proteins - just to rank transformed first 
clock <- names(d1)[c(160,161)] 
write.csv(d1, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/d1.csv", row.names = F)

####################################################################################

# Read in the prepped cases files for each trait to be assessed 

Diabetes <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/Diabetes_combined.csv")

Stroke <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/Stroke_combined.csv")

IHD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/IHD_combined.csv")

DEM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_Cox_proteins_20k/Trait_preps/all_dementia.csv")

####################################################################################

# Subset to waves

target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")

d1 <- left_join(d1, target, by = "Sample_Name")

W1 <- d1[which(d1$Set %in% "wave1"),]
W3 <- d1[which(d1$Set %in% "wave3"),]
W4 <- d1[which(d1$Set %in% "wave4"),]

d1 <- W4

####################################################################################

## HAZARD MODELS - PREDICTING TIME-TO-ONSET OF DISEASES FROM STUDY BASELINE (2006)

####################################################################################

d1_DEM <- d1

mat_hazard_ad <- matrix(nrow=length(clock),ncol=9)
output_hazard_DEM<- as.data.frame(mat_hazard_ad)
for(j in 1:length(clock)){ 
  tryCatch({ 
  dat1= d1_DEM
  tmp1 = DEM[which(DEM$id %in% dat1$Sample_Name),] # 135 individuals 
  
  ## Obtain Age of Onset 
  affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
  age_onset = DEM[,c("first", "id")]
  affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
  affected$Event = 1
  affected$yoe = substring(affected$first, 1, 4)
  affected$moe = substring(affected$first, 5,6)
  affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
  affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
  affected$age_event = affected$age_event1 + affected$month_event1
  affected$first = NULL
  affected$yoe = NULL 
  affected$moe = NULL
  affected$month_event1 = NULL 
  affected$age_event1 = NULL
  
  healthy = dat1[-which(dat1$Sample_Name %in% DEM$id),]
  healthy$Event = 0
  healthy$age_event = 0 
  affected$id.y <- NULL
  healthy$id <- NULL
  names(affected)[names(affected)=="id"] <- "Sample_Name"
  cox = rbind(affected, healthy)
  
  cox$age_death = 0
  cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
  cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
  cox$tte = cox$age_at_event - cox$age.x
  cox$tte = as.numeric(cox$tte)
  cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
  cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
  cox$Event = as.numeric(cox$Event)
  cox$tte<-as.numeric(cox$tte)

  cox <- cox[complete.cases(cox[,clock[[j]]]),]
  
  cox = cox[cox$age_at_event >=65,]
  mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + factor(cox$sex.x) + cox$age.x + (1|cox$Sample_Name), varlist = kin_model*2)
  print(clock[[j]])
  print("DEM")
  output_hazard_DEM[j,1] <- as.character(clock[[j]])
  output_hazard_DEM[j,2] <- as.character("Dementia")
  output_hazard_DEM[j,3:5]<-round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
  output_hazard_DEM[j,6] <- extract_coxme_table(mod)[1,4]
  output_hazard_DEM[j,7] <- mod$n[1]
  output_hazard_DEM[j,8] <- mod$n[2]-mod$n[1]
  output_hazard_DEM[j,9] <-cox.zph(mod)[1][[1]][3]

  # Get time to event info for included cases, with max tte reported
  cox$Event <- ifelse(cox$tte < 0, "NA", cox$Event) # add this in to change event to NA as well as tte
  check <- cox %>% filter(Event == "1")
  check <- check %>% filter(!tte == 'NA')

  mean_tte <- mean(check$tte, na.rm = T) %>% round(digits = 1)
  sd_tte <- sd(check$tte, na.rm = T) %>% round(digits = 1) 
  mean_sd <- paste0(mean_tte, " ", "(", sd_tte, ")") 
  output_hazard_DEM[j,10] <- mean_sd

  max <- max(check$tte, na.rm = T)
  output_hazard_DEM[j,11] <- max

  write.csv(cox, paste0("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/Cox_DEM_", clock[[j]], ".csv"), row.names = F)

  }, error = function(e) cat("skipped"))
} 


## Create List of Remaining Dataframes

my.list = list(Stroke,Diabetes,IHD)
names = list("Stroke","Diabetes","IHD")
names(my.list) <- names 

l=lapply(my.list, "[", c(1:2))

names(d1)[names(d1) == "heart_disease_Y"] <- "IHD"
names(d1)[names(d1) == "stroke_Y"] <- "Stroke"
names(d1)[names(d1) == "diabetes_Y"] <- "Diabetes"

mat_hazard <- matrix(nrow=200*length(my.list),ncol=9)
output_hazard <- as.data.frame(mat_hazard)
k=c(0,200,400,600,800,1000,1200,1400)

## Loop of Survival Models - Longitudinal Associations
for(j in 1:length(clock)){
  for(i in 1:length(l)){ 
   tryCatch({ 
    tmp <- l[[i]]
    
    ## Exclude Indiviudals who Reported Disease at Study Baseline  
    dat1= d1[-which(d1[,names[[i]]] %in% 1),]
    tmp1 = tmp[which(tmp$id %in% dat1$Sample_Name),]
    
    ## Obtain Age of Onset 
    affected = dat1[which(dat1$Sample_Name %in% tmp1$id),] 
    age_onset = tmp[,c("first", "id")]
    affected = merge(age_onset, affected, by.x = "id", by.y = "Sample_Name")
    affected$Event = 1
    affected$yoe = substring(affected$first, 1, 4)
    affected$moe = substring(affected$first, 5,6)
    affected$month_event1 = (as.numeric(affected$moe) - as.numeric(affected$mob))/12
    affected$age_event1 = as.numeric(affected$yoe) - as.numeric(affected$yob)
    affected$age_event = affected$age_event1 + affected$month_event1
    affected$first = NULL
    affected$yoe = NULL 
    affected$moe = NULL
    affected$month_event1 = NULL 
    affected$age_event1 = NULL
    
    healthy = dat1[-which(dat1$Sample_Name %in% tmp$id),]
    healthy$Event = 0
    healthy$age_event = 0 
    affected$id.y <- NULL
    healthy$id <- NULL
    names(affected)[names(affected)=="id"] <- "Sample_Name"
    cox = rbind(affected, healthy)
    
    cox$age_death = 0
    cox$age_death = ifelse(cox$dead %in% 1, cox$aged, 0)
    cox$age_at_event = ifelse(cox$Event %in% 1, cox$age_event, (ifelse(cox$dead %in% 1 & cox$Event %in% 0, cox$age_death, cox$aged)))
    cox$tte = cox$age_at_event - cox$age.x
    cox$tte = as.numeric(cox$tte)
    cox$tte <- ifelse(cox$tte < -1, "NA", cox$tte)
    cox$tte = ifelse(cox$tte < 0, 0, cox$tte)
    cox$Event = as.numeric(cox$Event)
    cox$tte<-as.numeric(cox$tte)
       
    cox <- cox[complete.cases(cox[,clock[[j]]]),] 

    mod = coxme(Surv(cox$tte, cox$Event) ~ scale(cox[,clock[[j]]]) + cox$age.x + factor(cox$sex.x) + (1|cox$Sample_Name), varlist = kin_model*2)
    print(names[[i]])
    print(clock[[j]])
    output_hazard[j+k[[i]],1] <- as.character(clock[[j]])
    output_hazard[j+k[[i]],2] <- as.character(names[[i]])
    output_hazard[j+k[[i]],3:5] <- round(exp(cbind(coef(mod), confint(mod)))[1,1:3],2)
    output_hazard[j+k[[i]],6] <- extract_coxme_table(mod)[1,4]
    output_hazard[j+k[[i]],7] <- mod$n[1]
    output_hazard[j+k[[i]],8] <- mod$n[2]-mod$n[1]
    output_hazard[j+k[[i]],9] <-cox.zph(mod)[1][[1]][3]

    cox$Event <- ifelse(cox$tte < 0, "NA", cox$Event) # add this in to change event to NA as well as tte
    check <- cox %>% filter(Event == "1")
    check <- check %>% filter(!tte == 'NA')

    mean_tte <- mean(check$tte, na.rm = T) %>% round(digits = 1)
    sd_tte <- sd(check$tte, na.rm = T) %>% round(digits = 1) 
    mean_sd <- paste0(mean_tte, " ", "(", sd_tte, ")") 
    output_hazard[j+k[[i]],10]  <- mean_sd

    max <- max(check$tte, na.rm = T)
    output_hazard[j+k[[i]],11] <- max

    write.csv(cox, paste0("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/Cox_", names[[i]], "_", clock[[j]], ".csv"), row.names = F)

   }, error = function(e) cat("skipped"))
  } 
} 



##################################################################################################

### Filtering by relatedness 

# Get cases for basic models in W4 and identify related individuals to remove from the training sample

library(tidyverse)

## GDF15 

# Get case and control lists for the basic models in W4 joint to ped file to source famids
ped <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree_formatted.csv")
names(ped)[2] <- "Sample_Name"

diab <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/Cox_Diabetes_gdf15.csv") # 8142
diab <- diab %>% filter(!tte == "NA")
diab <- left_join(diab, ped, by = "Sample_Name")
diab_cases <- diab[which(diab$Event == "1"),] 
diab_con <- diab[which(diab$Event == "0"),] 

IHD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/Cox_IHD_gdf15.csv")
IHD <- left_join(IHD, ped, by = "Sample_Name")
IHD <- IHD %>% filter(!tte == "NA") 
IHD_cases <- IHD[which(IHD$Event == "1"),] 
IHD_con <- IHD[which(IHD$Event == "0"),] 

stroke <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/Cox_Stroke_gdf15.csv")
stroke <- left_join(stroke, ped, by = "Sample_Name")
stroke <- stroke %>% filter(!tte == "NA") 
stroke_cases <- stroke[which(stroke$Event == "1"),] 
stroke_con <- stroke[which(stroke$Event == "0"),] 

DEM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/Cox_DEM_gdf15.csv")
DEM <- left_join(DEM, ped, by = "Sample_Name")
DEM <- DEM %>% filter(!tte == "NA") 
DEM_cases <- DEM[which(DEM$Event == "1"),]
DEM_con <- DEM[which(DEM$Event == "0"),]

# Now identify those who are related to cases for each disease in the W1/W3 training sample - i.e. to remove from training 
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
W1W3 <- target[which(target$Set %in% c("wave1", "wave3")),] # 9537
W1W3 <- left_join(W1W3, ped, by = "Sample_Name")

# Index individuals for exlcusion in training 
ex_diab <- W1W3[which(W1W3$famid %in% diab_cases$famid),]
ex_IHD <- W1W3[which(W1W3$famid %in% IHD_cases$famid),]
ex_stroke <- W1W3[which(W1W3$famid %in% stroke_cases$famid),]
ex_DEM <- W1W3[which(W1W3$famid %in% DEM_cases$famid),]

bind <- rbind(ex_diab, ex_IHD)
bind <- rbind(bind, ex_stroke)
bind <- rbind(bind, ex_DEM)

length(unique(bind$Sample_Name))

bind <- bind[c(1,2)]
write.csv(bind, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/GDF_exclusion_training_W1W3.csv", row.names = F)

# Check to see what the training set is after exclusions
W1W3_train <- W1W3[-which(W1W3$Sample_Name %in% bind$Sample_Name),]

# Check to make sure no relatives of cases are in the training set 
which(W1W3_train$famid %in% diab_cases$famid)
which(W1W3_train$famid %in% stroke_cases$famid)
which(W1W3_train$famid %in% IHD_cases$famid)
which(W1W3_train$famid %in% DEM_cases$famid)

## Now subset W4 to those unrelated in to W1W3 training set as test set 
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
W4 <- target[which(target$Set %in% c("wave4")),] # 8876
W4 <- left_join(W4, ped, by = "Sample_Name")

W4_test <- W4[-which(W4$famid %in% W1W3_train$famid),]

# Check the cases are still present in the test set 
dim(diab_cases)
dim(stroke_cases)
dim(IHD_cases)
dim(DEM_cases)

length(which(W4_test$Sample_Name %in% diab_cases$Sample_Name))
length(which(W4_test$Sample_Name %in% stroke_cases$Sample_Name))
length(which(W4_test$Sample_Name %in% IHD_cases$Sample_Name))
length(which(W4_test$Sample_Name %in% DEM_cases$Sample_Name))


# Check to see relatedness to cases in test set 
length(which(W4_test$famid %in% diab_cases$famid))
length(which(W4_test$famid %in% stroke_cases$famid))
length(which(W4_test$famid %in% IHD_cases$famid))
length(which(W4_test$famid %in% DEM_cases$famid))


# Write out these as test sets of W4 individuals unrelated to those in W1W3 training 

write.csv(W4_test, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/unrelated_W4_GDF.csv", row.names = F)


################################################################################################################

### BNP

library(tidyverse)

# Get case and control lists for the basic models in W4 joint to ped file to source famids
ped <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/Proxies_stradl/KORA_train_recieved_020221/Cox/pedigree_formatted.csv")
names(ped)[2] <- "Sample_Name"

diab <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/Cox_Diabetes_nt.probnp.csv") 
diab <- diab %>% filter(!tte == "NA") 
diab <- left_join(diab, ped, by = "Sample_Name")
diab_cases <- diab[which(diab$Event == "1"),] 
diab_con <- diab[which(diab$Event == "0"),] 

IHD <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/Cox_IHD_nt.probnp.csv")
IHD <- left_join(IHD, ped, by = "Sample_Name")
IHD <- IHD %>% filter(!tte == "NA") 
IHD_cases <- IHD[which(IHD$Event == "1"),]
IHD_con <- IHD[which(IHD$Event == "0"),]

stroke <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/Cox_Stroke_nt.probnp.csv")
stroke <- left_join(stroke, ped, by = "Sample_Name")
stroke <- stroke %>% filter(!tte == "NA") 
stroke_cases <- stroke[which(stroke$Event == "1"),]
stroke_con <- stroke[which(stroke$Event == "0"),]

DEM <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/Cox_DEM_nt.probnp.csv")
DEM <- left_join(DEM, ped, by = "Sample_Name")
DEM <- DEM %>% filter(!tte == "NA") 
DEM_cases <- DEM[which(DEM$Event == "1"),] 
DEM_con <- DEM[which(DEM$Event == "0"),]


# Now identify those who are related to cases for each disease in the W1/W3 training sample - i.e. to remove from training 
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
W1W3 <- target[which(target$Set %in% c("wave1", "wave3")),] # 9537
W1W3 <- left_join(W1W3, ped, by = "Sample_Name")

# Index individuals for exlcusion in training 
ex_diab <- W1W3[which(W1W3$famid %in% diab_cases$famid),]
ex_IHD <- W1W3[which(W1W3$famid %in% IHD_cases$famid),]
ex_stroke <- W1W3[which(W1W3$famid %in% stroke_cases$famid),]
ex_DEM <- W1W3[which(W1W3$famid %in% DEM_cases$famid),]

bind <- rbind(ex_diab, ex_IHD)
bind <- rbind(bind, ex_stroke)
bind <- rbind(bind, ex_DEM)

bind <- bind[c(1,2)]
write.csv(bind, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/BNP_exclusion_training_W1W3.csv", row.names = F)


# Check to see what the training set is after exclusions
W1W3_train <- W1W3[-which(W1W3$Sample_Name %in% bind$Sample_Name),]

# Check to make sure no relatives of cases are in the training set 
which(W1W3_train$famid %in% diab_cases$famid)
which(W1W3_train$famid %in% stroke_cases$famid)
which(W1W3_train$famid %in% IHD_cases$famid)
which(W1W3_train$famid %in% DEM_cases$famid)

## Now subset W4 to those unrelated in to W1W3 training set as test set 
target <- readRDS("/Cluster_Filespace/Marioni_Group/GS/GS_methylation/GS20k/GS20k_Targets.rds")
W4 <- target[which(target$Set %in% c("wave4")),] # 8876
W4 <- left_join(W4, ped, by = "Sample_Name")

W4_test <- W4[-which(W4$famid %in% W1W3_train$famid),]

# Check the cases are still present in the test set 
dim(diab_cases)
dim(stroke_cases)
dim(IHD_cases)
dim(DEM_cases)

length(which(W4_test$Sample_Name %in% diab_cases$Sample_Name))
length(which(W4_test$Sample_Name %in% stroke_cases$Sample_Name))
length(which(W4_test$Sample_Name %in% IHD_cases$Sample_Name))
length(which(W4_test$Sample_Name %in% DEM_cases$Sample_Name))

# Check to see relatedness to cases in test set 
length(which(W4_test$famid %in% diab_cases$famid))
length(which(W4_test$famid %in% stroke_cases$famid))
length(which(W4_test$famid %in% IHD_cases$famid))
length(which(W4_test$famid %in% DEM_cases$famid))

# Write out these as test sets of W4 individuals unrelated to those in W1W3 training 

write.csv(W4_test, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/unrelated_W4_BNP.csv", row.names = F)
