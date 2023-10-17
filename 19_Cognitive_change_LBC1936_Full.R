#### library ####
library("tidyverse")
library("lavaan")
library("readr")

#### data #####

scores <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/LBC_episcores_adjusted_technical.csv")
LBC <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/LBC1936_HannahSmith.csv")
covars <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/LBC_processed_data.csv")

#select columns 
data <- LBC %>% select(lbc36no, lnseq_w1:wtar_total_w4,lnseq_w5:wtar_total_w5)


#==================================================================
#     Select variables and recode missing
#==================================================================

#recode missing
data[data == -999] <- NA
data[data == -777] <- NA
data[data == 999] <- NA



#check_symsear_w1 <- data %>% filter(symsear_w1 < 0) %>% select(symsear_w1) # one person has -1, one person has -4 
#check_symsear_w4 <- data %>% filter(symsear_w4 <0) %>% select(symsear_w4) # one person has -3 

data$symsear_w1 =ifelse(data$symsear_w1 < 0, NA, data$symsear_w1)
data$symsear_w4 =ifelse(data$symsear_w4 %in% -3, NA, data$symsear_w4)

#### join data #####
data <- data %>% left_join(covars, by = "lbc36no")
data <- data %>% left_join(scores, by = "Basename")


dset_mod <- mutate(data,
                   blkdes_w1 = blkdes_w1/2,
                   blkdes_w2 = blkdes_w2/2,
                   blkdes_w3 = blkdes_w3/2,
                   blkdes_w4 = blkdes_w4/2,
                   blkdes_w5 = blkdes_w5/2,
                   vftot_w1 = vftot_w1/2,
                   vftot_w2 = vftot_w2/2,
                   vftot_w3 = vftot_w3/2,
                   vftot_w4 = vftot_w4/2,
                   vftot_w5 = vftot_w5/2,
                   lmtotal_w1 = lmtotal_w1/3,
                   lmtotal_w2 = lmtotal_w2/3,
                   lmtotal_w3 = lmtotal_w3/3,
                   lmtotal_w4 = lmtotal_w4/3,
                   lmtotal_w5 = lmtotal_w5/3,
                   digback_w1 = 3*digback_w1,
                   digback_w2 = 3*digback_w2,
                   digback_w3 = 3*digback_w3,
                   digback_w4 = 3*digback_w4,
                   digback_w5 = 3*digback_w5,
                   digsym_w1 = digsym_w1/2,
                   digsym_w2 = digsym_w2/2,
                   digsym_w3 = digsym_w3/2,
                   digsym_w4 = digsym_w4/2,
                   digsym_w5 = digsym_w5/2,
                   ittotal_w1 = ittotal_w1/2,
                   ittotal_w2 = ittotal_w2/2,
                   ittotal_w3 = ittotal_w3/2,
                   ittotal_w4 = ittotal_w4/2,
                   ittotal_w5 = ittotal_w5/2,
                   crtmean_w1 = -50 * crtmean_w1,
                   crtmean_w2 = -50 * crtmean_w2,
                   crtmean_w3 = -50 * crtmean_w3,
                   crtmean_w4 = -50 * crtmean_w4,
                   crtmean_w5 = -50 * crtmean_w5, 
                   depind_w1 = depind_w1/1000)

dset_mod <- as.data.frame(dset_mod)
dset_mod$sex <- as.factor(dset_mod$sex)
dset_mod$smokcat_w1 <- as.factor(dset_mod$smokcat_w1)

#### growth curve models for individual cognitive tests ####

pgmodel1 <- '
Imatreas =~ 1*matreas_w1 + 1*matreas_w2 + 1*matreas_w3 + 1*matreas_w4 + 1*matreas_w5
Smatreas =~ 0*matreas_w1 + 2.98*matreas_w2 + 6.75*matreas_w3 + 9.82*matreas_w4 + 12.54*matreas_w5
'
fit1 <- growth(pgmodel1, dset_mod, missing = "ml.x")
summary(fit1, standardized = T)

pgmodel2 <-'
Iblkdes =~ 1*blkdes_w1 + 1*blkdes_w2 + 1*blkdes_w3 + 1*blkdes_w4 + 1*blkdes_w5
Sblkdes=~ 0*blkdes_w1 + 2.98*blkdes_w2 + 6.75*blkdes_w3 + 9.82*blkdes_w4 + 12.54*blkdes_w5
'
fit2 <- growth(pgmodel2, dset_mod, missing = "ml.x")
summary(fit2, standardized = T)

pgmodel3 <- '
Ispantot =~ 1*spantot_w1 + 1*spantot_w2 + 1*spantot_w3 + 1*spantot_w4 + 1*spantot_w5
Sspantot=~ 0*spantot_w1 + 2.98*spantot_w2 + 6.75*spantot_w3 + 9.82*spantot_w4 + 12.54*spantot_w5
'
fit3 <- growth(pgmodel3, dset_mod, missing = "ml.x")
summary(fit3, standardized = T)

pgmodel4 <- '
Inart =~ 1*nart_w1 + 1*nart_w2 + 1*nart_total_w3 + 1*nart_total_w4 + 1*nart_total_w5
Snart =~ 0*nart_w1 + 2.98*nart_w2 + 6.75*nart_total_w3 + 9.82*nart_total_w4 + 12.54*nart_total_w5
'
fit4 <- growth(pgmodel4, dset_mod, missing = "ml.x")
summary(fit4, standardized = T)

pgmodel5 <- '
Iwtar =~ 1*wtar_w1 + 1*wtar_w2 + 1*wtar_total_w3 + 1*wtar_total_w4 + 1*wtar_total_w5
Swtar =~ 0*wtar_w1 + 2.98*wtar_w2 + 6.75*wtar_total_w3 + 9.82*wtar_total_w4 + 12.54*wtar_total_w5
'
fit5 <- growth(pgmodel5, dset_mod, missing = "ml.x")
summary(fit5, standardized = T)

pgmodel6 <- '
Ivftot =~ 1*vftot_w1 + 1*vftot_w2 + 1*vftot_w3 + 1*vftot_w4 + 1*vftot_w5
Svftot =~ 0*vftot_w1 + 2.98*vftot_w2 + 6.75*vftot_w3 + 9.82*vftot_w4 + 12.54*vftot_w5
'
fit6 <- growth(pgmodel6, dset_mod, missing = "ml.x")
summary(fit6, standardized = T)

pgmodel7 <-'
Ivpatotal =~ 1*vpatotal_w1 + 1*vpatotal_w2 + 1*vpatotal_w3 + 1*vpatotal_w4 + 1*vpa_total_w5
Svpatotal =~ 0*vpatotal_w1 + 2.98*vpatotal_w2 + 6.75*vpatotal_w3 + 9.82*vpatotal_w4 + 12.54*vpa_total_w5
'
fit7 <- growth(pgmodel7, dset_mod, missing = "ml.x")
summary(fit7, standardized = T)

pgmodel8 <- '
Ilmtotal =~ 1*lmtotal_w1 + 1*lmtotal_w2 + 1*lmtotal_w3 + 1*lmtotal_w4 + 1*lmtotal_w5
Slmtotal =~ 0*lmtotal_w1 + 2.98*lmtotal_w2 + 6.75*lmtotal_w3 + 9.82*lmtotal_w4 + 12.54*lmtotal_w5
'
fit8 <- growth(pgmodel8, dset_mod, missing = "ml.x")
summary(fit8, standardized = T)

pgmodel9 <- '
Idigback =~ 1*digback_w1 + 1*digback_w2 + 1*digback_w3 + 1*digback_w4 + 1*digback_w5
Sdigback =~ 0*digback_w1 + 2.98*digback_w2 + 6.75*digback_w3 + 9.82*digback_w4 + 12.54*digback_w5
'
fit9 <- growth(pgmodel9, dset_mod, missing = "ml.x")
summary(fit9, standardized = T)

pgmodel10 <- '
Isymsear =~ 1*symsear_w1 + 1*symsear_w2 + 1*symsear_w3 + 1*symsear_w4 + 1*symsear_w5
Ssymsear =~ 0*symsear_w1 + 2.98*symsear_w2 + 6.75*symsear_w3 + 9.82*symsear_w4 + 12.54*symsear_w5
'
fit10 <- growth(pgmodel10, dset_mod, missing = "ml.x")
summary(fit10, standardized = T)

pgmodel11 <- '
Idigsym =~ 1*digsym_w1 + 1*digsym_w2 + 1*digsym_w3 + 1*digsym_w4 + 1*digsym_w5
Sdigsym =~ 0*digsym_w1 + 2.98*digsym_w2 + 6.75*digsym_w3 + 9.82*digsym_w4 + 12.54*digsym_w5
'
fit11 <- growth(pgmodel11, dset_mod, missing = "ml.x")
summary(fit11, standardized = T)

pgmodel12 <- '
Iittotal =~ 1*ittotal_w1 + 1*ittotal_w2 + 1*ittotal_w3 + 1*ittotal_w4 + 1*ittotal_w5
Sittotal =~ 0*ittotal_w1 + 2.98*ittotal_w2 + 6.75*ittotal_w3 + 9.82*ittotal_w4 + 12.54*ittotal_w5
'
fit12 <- growth(pgmodel12, dset_mod, missing = "ml.x")
summary(fit12, standardized = T)

pgmodel13 <- '
Icrtmean =~ 1*crtmean_w1 + 1*crtmean_w2 + 1*crtmean_w3 + 1*crtmean_w4 + 1*crtmean_w5
Scrtmean =~ 0*crtmean_w1 + 2.98*crtmean_w2 + 6.75*crtmean_w3 + 9.82*crtmean_w4 + 12.54*crtmean_w5
'
fit13 <- growth(pgmodel13, dset_mod, missing = "ml.x")
summary(fit13, standardized = T)



general_4p <- '

#latent variables 
Ig =~  Iblkdes + Imatreas  + Ispantot + Inart + Iwtar + Ivftot + Ivpatotal + Ilmtotal +
  Idigback + Isymsear + Idigsym + Iittotal + Icrtmean

Sg =~ Sblkdes + Smatreas + Sspantot + Snart + Swtar + Svftot + Svpatotal + Slmtotal +
  Sdigback + Ssymsear + Sdigsym + Sittotal + Scrtmean


#indicator as scaling reference: loading=1, int=0
Iblkdes ~ 0*1
Sblkdes ~ 0*1 

#fixed
Sspantot ~~ 0*Sspantot
Sdigback ~~ 0*Sdigback
Smatreas ~~ 0*Smatreas


#covariances 
Iblkdes ~~ Imatreas 
Iblkdes ~~ Ispantot
Imatreas ~~ Ispantot

#Sblkdes ~~ Smatreas 
#Sblkdes ~~ Sspantot
#Smatreas ~~ Sspantot

Iwtar ~~ Inart 
Iwtar ~~ Ivftot
Inart ~~ Ivftot

Swtar ~~ Snart 
Swtar ~~ Svftot
Snart ~~ Svftot

Ilmtotal ~~ Ivpatotal
Ilmtotal ~~ Idigback
Ivpatotal ~~ Idigback

Slmtotal ~~ Svpatotal
#Slmtotal ~~ Sdigback
#Svpatotal ~~ Sdigback

Iittotal ~~ Idigsym 
Iittotal ~~ Isymsear
Iittotal ~~ Icrtmean
Idigsym ~~ Isymsear
Idigsym ~~ Icrtmean
Isymsear ~~ Icrtmean

Sittotal ~~ Sdigsym 
Sittotal ~~ Ssymsear
Sittotal ~~ Scrtmean
Sdigsym ~~ Ssymsear
Sdigsym ~~ Scrtmean
Ssymsear ~~ Scrtmean
'

fitGen_4p <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, pgmodel5, pgmodel6, pgmodel7, pgmodel8, pgmodel9, pgmodel10, pgmodel11, pgmodel12, pgmodel13, general_4p), dset_mod,  missing = "ml.x")
fitmeasures(fitGen_4p, c("cfi", "tli", "RMSEA", "SRMR")) # cfi   tli rmsea  srmr : 0.958 0.957 0.029 0.062
summary(fitGen_4p, standardized = T)

list_e <- colnames(dset_mod[103:104])

results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)

rownames(results) = results$SeqId # This will allow you to index with results[i,]
episcores = colnames(dset_mod)[which(colnames(dset_mod)%in% list_e)]
for(i in episcores) { 
  dset_mod$tmp = unlist(dset_mod[,i])
  
  reg_Ig <- '
 Ig ~ tmp + ageyearsw1 + sex + neut + lymph + mono + eosin + baso + smokingScore + bmi_w1 + depind_w1 + alcunitwk_w1
'
  
  fitreg_g <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, pgmodel5, pgmodel6, pgmodel7, pgmodel8, pgmodel9, pgmodel10, pgmodel11, pgmodel12, pgmodel13,general_4p, reg_Ig), dset_mod,  missing = "ml.x")
  
  output = standardizedSolution(fitreg_g)
  ind = which(output$lhs == "Ig" & output$op == "~" & output$rhs=="tmp")
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(fitreg_g)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(fitreg_g, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  phenotype <- output[ind, "lhs"]
  
  
  results[i,1] <- i
  results[i,2] <- n
  results[i,3] <- Beta
  results[i,4] <- SE
  results[i,5] <- p
  results[i,6] <- ci.upper
  results[i,7] <- ci.lower
  results[i,8] <- cfi
  results[i,9] <- rmsea
  results[i,10] <- srmr
  results[i,11] <- tli
  results[i,12] <- phenotype
  
}
write.csv(results, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/General_Cog_ability_Score_intercept_LBC1936_assocs_2_EpiScores_full_15122022.csv", row.names = FALSE)


list_e <- colnames(dset_mod[103:104])

results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)

rownames(results) = results$SeqId # This will allow you to index with results[i,]
episcores = colnames(dset_mod)[which(colnames(dset_mod)%in% list_e)]
for(i in episcores) { 
  dset_mod$tmp = unlist(dset_mod[,i])
  
  reg_Sg <- '
 Sg ~ tmp + ageyearsw1 + sex + neut + lymph + mono + eosin + baso + smokingScore + bmi_w1 + depind_w1 + alcunitwk_w1
'
  
  fitreg_g <- growth(model = c(pgmodel1, pgmodel2, pgmodel3, pgmodel4, pgmodel5, pgmodel6, pgmodel7, pgmodel8, pgmodel9, pgmodel10, pgmodel11, pgmodel12, pgmodel13,general_4p, reg_Sg), dset_mod,  missing = "ml.x")
  
  output = standardizedSolution(fitreg_g)
  ind = which(output$lhs == "Sg" & output$op == "~" & output$rhs=="tmp")
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(fitreg_g)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(fitreg_g, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  phenotype <- output[ind, "lhs"]
  
  
  results[i,1] <- i
  results[i,2] <- n
  results[i,3] <- Beta
  results[i,4] <- SE
  results[i,5] <- p
  results[i,6] <- ci.upper
  results[i,7] <- ci.lower
  results[i,8] <- cfi
  results[i,9] <- rmsea
  results[i,10] <- srmr
  results[i,11] <- tli
  results[i,12] <- phenotype
  
}

write.csv(results, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/General_Cog_change(slope)_LBC1936_assocs_2_EpiScores_full_15122022.csv", row.names = FALSE)
