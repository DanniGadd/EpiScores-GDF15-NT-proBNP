library(readr)
library(lavaan)
library(tidyverse)

scores <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/20k_GS_training/LBC_episcores_adjusted_technical.csv")
LBC <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/LBC1936_HannahSmith.csv")
covars <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/LBC_processed_data.csv")

## Create age in years
LBC$Age_w1 <- LBC$agedays_w1 / 365.25
LBC$Age_w2 <- LBC$agedays_w2 / 365.25
LBC$Age_w3 <- LBC$agedays_w3 / 365.25
LBC$Age_w4 <- LBC$agedays_w4 / 365.25
LBC$Age_w5 <- LBC$agedays_w5 / 365.25

attach(LBC)
LBC$wmh_mm3_log_w2 <- log(wmh_mm3_w2 + 1)
LBC$wmh_mm3_log_w3 <- log(wmh_mm3_w3 + 1)
LBC$wmh_mm3_log_w4 <- log(wmh_mm3_w4 + 1)
LBC$wmh_mm3_log_w5 <- log(wmh_mm3_w5 + 1)
detach()

data <- LBC %>% select(lbc36no, brain_mm3_w2:brain_mm3_w5, gm_mm3_w2:nawm_mm3_w5, Age_w1:wmh_mm3_log_w5, ICV_mm3_wX)


data[data == -999] <- NA
data[data == -777] <- NA
data[data == 999] <- NA

#### join data #####
data <- data %>% left_join(covars, by = "lbc36no")
data <- data %>% left_join(scores, by = "Basename")

#
data$sex <- as.factor(data$sex)



data <- mutate(data,
               ICV_mm3_wX = scale(ICV_mm3_wX),
               Age_w1 = scale(Age_w1),
               Age_w2 = scale(Age_w2),
               Age_w3 = scale(Age_w3),
               Age_w4 = scale(Age_w4),
               Age_w5 = scale(Age_w5),
               depind_w1 = scale(depind_w1))




data$gm_mm3_w2 <- data$gm_mm3_w2/10000
data$gm_mm3_w3 <- data$gm_mm3_w3/10000
data$gm_mm3_w4 <- data$gm_mm3_w4/10000
data$gm_mm3_w5 <- data$gm_mm3_w5/10000

data$brain_mm3_w2 <- data$brain_mm3_w2/10000
data$brain_mm3_w3 <- data$brain_mm3_w3/10000
data$brain_mm3_w4 <- data$brain_mm3_w4/10000
data$brain_mm3_w5 <- data$brain_mm3_w5/10000

data$nawm_mm3_w2 <- data$nawm_mm3_w2/10000
data$nawm_mm3_w3 <- data$nawm_mm3_w3/10000
data$nawm_mm3_w4 <- data$nawm_mm3_w4/10000
data$nawm_mm3_w5 <- data$nawm_mm3_w5/10000



########################################################################################
########################################################################################

GM_model <- '
   # latent variables
     IGMV =~ 1*gm_mm3_w2 + 1*gm_mm3_w3 + 1*gm_mm3_w4 + 1*gm_mm3_w5
     SGMV =~ 0*gm_mm3_w2 + 3.76*gm_mm3_w3 + 6.83*gm_mm3_w4 + 9.55*gm_mm3_w5
  
   # regressions 
     gm_mm3_w2 ~ Age_w2
     gm_mm3_w3 ~ Age_w3
     gm_mm3_w4 ~ Age_w4
     gm_mm3_w5 ~ Age_w5
    
'
fit_GM <- growth(model = GM_model, data  = data, missing = "ml.x")
summary(fit_GM, standardized = TRUE, fit.measures = TRUE)
fitmeasures(fit_GM, c("cfi", "tli", "RMSEA", "SRMR"))



list_e <- colnames(data)[55:56]

results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)
results_slope <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)

rownames(results) = results$SeqId # This will allow you to index with results[i,]
rownames(results_slope) = results_slope$SeqId

episcores = colnames(data)[which(colnames(data)%in% list_e)]
for(i in episcores) { 
  data$tmp = unlist(data[,i])
  
  
  GM_model_reg <-  '
    # latent variables
    IGMV =~ 1*gm_mm3_w2 + 1*gm_mm3_w3 + 1*gm_mm3_w4 + 1*gm_mm3_w5
    SGMV =~ 0*gm_mm3_w2 + 3.76*gm_mm3_w3 + 6.83*gm_mm3_w4 + 9.55*gm_mm3_w5
    
    # regressions 
    gm_mm3_w2 ~ Age_w2
    gm_mm3_w3 ~ Age_w3
    gm_mm3_w4 ~ Age_w4
    gm_mm3_w5 ~ Age_w5
    tmp ~ Age_w1
    IGMV ~ tmp  + ICV_mm3_wX + sex + smokingScore + bmi_w1 + depind_w1 + alcunitwk_w1
    SGMV ~ tmp + sex + smokingScore + bmi_w1 + depind_w1 + alcunitwk_w1
    
    '
  
  GM_fit_reg <- growth(model = GM_model_reg, data = data, missing = "ml.x" )
  
  
  output = standardizedSolution(GM_fit_reg)
  ind = which(output$lhs == "IGMV" & output$op == "~" & output$rhs=="tmp")
  snd = which(output$lhs == "SGMV" & output$op == "~" & output$rhs=="tmp")
  
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(GM_fit_reg)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(GM_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  pheno <- output[ind, "lhs"]
  
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
  results[i,12] <- pheno
  
  Beta_slope= output[snd, "est.std"] 
  SE_slope <- output[snd, "se"]
  p_slope <- output[snd, "pvalue"]
  n_slope <- nobs(GM_fit_reg)
  ci.upper_slope <- output[snd, "ci.upper"]
  ci.lower_slope <- output[snd, "ci.lower"]
  fitmeasures_slope <- fitMeasures(GM_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi_slope <- fitmeasures_slope["cfi"]
  rmsea_slope <- fitmeasures_slope["rmsea"]
  srmr_slope <- fitmeasures_slope["srmr"]
  tli_slope <- fitmeasures_slope["tli"]
  pheno_slope <- output[snd, "lhs"]
  
  results_slope[i,1] <- i
  results_slope[i,2] <- n_slope
  results_slope[i,3] <- Beta_slope
  results_slope[i,4] <- SE_slope
  results_slope[i,5] <- p_slope
  results_slope[i,6] <- ci.upper_slope
  results_slope[i,7] <- ci.lower_slope
  results_slope[i,8] <- cfi_slope
  results_slope[i,9] <- rmsea_slope
  results_slope[i,10] <- srmr_slope
  results_slope[i,11] <- tli_slope
  results_slope[i, 12] <- pheno_slope
  
  print(i)
  
  
}


write.csv(results, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/GM_intercept_LBC1936_assocs_2_EpiScores_mid_16122022.csv")
write.csv(results_slope, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/GM_slope_LBC1936_assocs_84_EpiScores_mid_16122022.csv")



###########################################################################################
###########################################################################################
Brain_model <- '
   # latent variables
     ITMV =~ 1*brain_mm3_w2 + 1*brain_mm3_w3 + 1*brain_mm3_w4 + 1*brain_mm3_w5
     STMV =~ 0*brain_mm3_w2 + 3.76*brain_mm3_w3 + 6.83*brain_mm3_w4 + 9.55*brain_mm3_w5
  
   # regressions 
     brain_mm3_w2 ~ Age_w2
     brain_mm3_w3 ~ Age_w3
     brain_mm3_w4 ~ Age_w4
     brain_mm3_w5 ~ Age_w5
    
'
Brain_fit <- growth(model = Brain_model, data  = data, missing = "ml.x")
summary(Brain_fit, standardized = TRUE, fit.measures = TRUE)
fitmeasures(Brain_fit, c("cfi", "tli", "RMSEA", "SRMR"))



list_e <- colnames(data)[55:56]

results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)
results_slope <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)

rownames(results) = results$SeqId # This will allow you to index with results[i,]
rownames(results_slope) = results_slope$SeqId

episcores = colnames(data)[which(colnames(data)%in% list_e)]
for(i in episcores) { 
  data$tmp = unlist(data[,i])
  
  
  Brain_model_reg <-  '
    # latent variables
     ITMV =~ 1*brain_mm3_w2 + 1*brain_mm3_w3 + 1*brain_mm3_w4 + 1*brain_mm3_w5
     STMV =~ 0*brain_mm3_w2 + 3.76*brain_mm3_w3 + 6.83*brain_mm3_w4 + 9.55*brain_mm3_w5
  
   # regressions 
     brain_mm3_w2 ~ Age_w2
     brain_mm3_w3 ~ Age_w3
     brain_mm3_w4 ~ Age_w4
     brain_mm3_w5 ~ Age_w5
     tmp ~ Age_w1
     ITMV ~ tmp + ICV_mm3_wX + sex + smokingScore + bmi_w1 + depind_w1 + alcunitwk_w1
     STMV ~ tmp + sex + smokingScore + bmi_w1 + depind_w1 + alcunitwk_w1
    '
  
  Brain_fit_reg <- growth(model = Brain_model_reg, data = data, missing = "ml.x" )
  
  
  output = standardizedSolution(Brain_fit_reg)
  ind = which(output$lhs == "ITMV" & output$op == "~" & output$rhs=="tmp")
  snd = which(output$lhs == "STMV" & output$op == "~" & output$rhs=="tmp")
  
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(Brain_fit_reg)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(Brain_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  pheno <- output[ind, "lhs"]
  
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
  results[i,12] <- pheno
  
  Beta_slope= output[snd, "est.std"] 
  SE_slope <- output[snd, "se"]
  p_slope <- output[snd, "pvalue"]
  n_slope <- nobs(Brain_fit_reg)
  ci.upper_slope <- output[snd, "ci.upper"]
  ci.lower_slope <- output[snd, "ci.lower"]
  fitmeasures_slope <- fitMeasures(Brain_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi_slope <- fitmeasures_slope["cfi"]
  rmsea_slope <- fitmeasures_slope["rmsea"]
  srmr_slope <- fitmeasures_slope["srmr"]
  tli_slope <- fitmeasures_slope["tli"]#
  pheno_slope <- output[snd, "lhs"]
  
  results_slope[i,1] <- i
  results_slope[i,2] <- n_slope
  results_slope[i,3] <- Beta_slope
  results_slope[i,4] <- SE_slope
  results_slope[i,5] <- p_slope
  results_slope[i,6] <- ci.upper_slope
  results_slope[i,7] <- ci.lower_slope
  results_slope[i,8] <- cfi_slope
  results_slope[i,9] <- rmsea_slope
  results_slope[i,10] <- srmr_slope
  results_slope[i,11] <- tli_slope
  results_slope[i,12] <- pheno_slope
  
  print(i)
  
}

write.csv(results, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/Total_Brain_intercept_LBC1936_assocs_2_EpiScores_mid_16122022.csv")
write.csv(results_slope, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/Total_Brain_slope_LBC1936_assocs_2_EpiScores_mid_16122022.csv")

###########################################################################################
###########################################################################################

NAWM_model <- '
   # latent variables
     INAWMV =~ 1*nawm_mm3_w2 + 1*nawm_mm3_w3 + 1*nawm_mm3_w4 + 1*nawm_mm3_w5
     SNAWMV =~ 0*nawm_mm3_w2 + 3.76*nawm_mm3_w3 + 6.83*nawm_mm3_w4 + 9.55*nawm_mm3_w5
  
   # regressions 
     nawm_mm3_w2 ~ Age_w2
     nawm_mm3_w3 ~ Age_w3
     nawm_mm3_w4 ~ Age_w4
     nawm_mm3_w5 ~ Age_w5
    
'
fit_NAWM <- growth(model = NAWM_model, data  = data, missing = "ml.x")
summary(fit_NAWM, standardized = TRUE, fit.measures = TRUE)
fitmeasures(fit_NAWM, c("cfi", "tli", "RMSEA", "SRMR"))



list_e <- colnames(data)[55:56]

results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)
results_slope <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)

rownames(results) = results$SeqId # This will allow you to index with results[i,]
rownames(results_slope) = results_slope$SeqId

episcores = colnames(data)[which(colnames(data)%in% list_e)]
for(i in episcores) { 
  data$tmp = unlist(data[,i])
  
  
  NAWM_model_reg <-  '
  # latent variables
  INAWMV =~ 1*nawm_mm3_w2 + 1*nawm_mm3_w3 + 1*nawm_mm3_w4 + 1*nawm_mm3_w5
  SNAWMV =~ 0*nawm_mm3_w2 + 3.76*nawm_mm3_w3 + 6.83*nawm_mm3_w4 + 9.55*nawm_mm3_w5
  
  # regressions 
  nawm_mm3_w2 ~ Age_w2
  nawm_mm3_w3 ~ Age_w3
  nawm_mm3_w4 ~ Age_w4
  nawm_mm3_w5 ~ Age_w5
  tmp ~ Age_w1
  INAWMV ~ tmp + ICV_mm3_wX + sex + smokingScore + bmi_w1 + depind_w1 + alcunitwk_w1
  SNAWMV ~ tmp + sex + smokingScore + bmi_w1 + depind_w1 + alcunitwk_w1
  '
  
  NAWM_fit_reg <- growth(model = NAWM_model_reg, data = data, missing = "ml.x" )
  
  
  output = standardizedSolution(NAWM_fit_reg)
  ind = which(output$lhs == "INAWMV" & output$op == "~" & output$rhs=="tmp")
  snd = which(output$lhs == "SNAWMV" & output$op == "~" & output$rhs=="tmp")
  
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(NAWM_fit_reg)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(NAWM_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  pheno <- output[ind, "lhs"]
  
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
  results[i,12] <- pheno
  
  Beta_slope= output[snd, "est.std"] 
  SE_slope <- output[snd, "se"]
  p_slope <- output[snd, "pvalue"]
  n_slope <- nobs(NAWM_fit_reg)
  ci.upper_slope <- output[snd, "ci.upper"]
  ci.lower_slope <- output[snd, "ci.lower"]
  fitmeasures_slope <- fitMeasures(NAWM_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi_slope <- fitmeasures_slope["cfi"]
  rmsea_slope <- fitmeasures_slope["rmsea"]
  srmr_slope <- fitmeasures_slope["srmr"]
  tli_slope <- fitmeasures_slope["tli"]
  pheno_slope <- output[snd, "lhs"]
  
  results_slope[i,1] <- i
  results_slope[i,2] <- n_slope
  results_slope[i,3] <- Beta_slope
  results_slope[i,4] <- SE_slope
  results_slope[i,5] <- p_slope
  results_slope[i,6] <- ci.upper_slope
  results_slope[i,7] <- ci.lower_slope
  results_slope[i,8] <- cfi_slope
  results_slope[i,9] <- rmsea_slope
  results_slope[i,10] <- srmr_slope
  results_slope[i,11] <- tli_slope
  results_slope[i,12] <- pheno_slope
  
  print(i)
}

write.csv(results, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/NAWM_intercept_LBC1936_assocs_2_EpiScores_mid_16122022.csv")
write.csv(results_slope, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/NAWM_slope_LBC1936_assocs_2_EpiScores_mid_16122022.csv")

###########################################################################################
###########################################################################################

WMH_model <- '
   # latent variables
     IWMHV =~ 1*wmh_mm3_log_w2 + 1*wmh_mm3_log_w3 + 1*wmh_mm3_log_w4 + 1*wmh_mm3_log_w5
     SWMHV =~ 0*wmh_mm3_log_w2 + 3.76*wmh_mm3_log_w3 + 6.83*wmh_mm3_log_w4 + 9.55*wmh_mm3_log_w5
  
   # regressions 
     wmh_mm3_log_w2 ~ Age_w2
     wmh_mm3_log_w3 ~ Age_w3
     wmh_mm3_log_w4 ~ Age_w4
     wmh_mm3_log_w5 ~ Age_w5
    
'
WMH_fit <- growth(model = WMH_model, data  = data, missing = "ml.x")
summary(WMH_fit, standardized = TRUE, fit.measures = TRUE)
fitmeasures(WMH_fit, c("cfi", "tli", "RMSEA", "SRMR"))



list_e <- colnames(data)[55:56]

results <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)
results_slope <- data.frame(SeqId = list_e, n = NA, beta = NA , SE = NA, P = NA, ci.upper = NA, ci.lower = NA, cfi = NA, rmsea = NA, srmr = NA, tli = NA, phenotype = NA)

rownames(results) = results$SeqId # This will allow you to index with results[i,]
rownames(results_slope) = results_slope$SeqId

episcores = colnames(data)[which(colnames(data)%in% list_e)]
for(i in episcores) { 
  data$tmp = unlist(data[,i])
  
  
  WMH_model_reg <-  '
  # latent variables
     IWMHV =~ 1*wmh_mm3_log_w2 + 1*wmh_mm3_log_w3 + 1*wmh_mm3_log_w4 + 1*wmh_mm3_log_w5
     SWMHV =~ 0*wmh_mm3_log_w2 + 3.76*wmh_mm3_log_w3 + 6.83*wmh_mm3_log_w4 + 9.55*wmh_mm3_log_w5
  
   # regressions 
     wmh_mm3_log_w2 ~ Age_w2
     wmh_mm3_log_w3 ~ Age_w3
     wmh_mm3_log_w4 ~ Age_w4
     wmh_mm3_log_w5 ~ Age_w5
     tmp ~ Age_w1
     IWMHV ~ tmp + ICV_mm3_wX + sex + smokingScore + bmi_w1 + depind_w1 + alcunitwk_w1
     SWMHV ~ tmp + sex + smokingScore + bmi_w1 + depind_w1 + alcunitwk_w1
  '
  
  WMH_fit_reg <- growth(model = WMH_model_reg, data = data, missing = "ml.x" )
  
  
  output = standardizedSolution(WMH_fit_reg)
  ind = which(output$lhs == "IWMHV" & output$op == "~" & output$rhs=="tmp")
  snd = which(output$lhs == "SWMHV" & output$op == "~" & output$rhs=="tmp")
  
  Beta = output[ind, "est.std"] 
  SE <- output[ind, "se"]
  p <- output[ind, "pvalue"]
  n <- nobs(WMH_fit_reg)
  ci.upper <- output[ind, "ci.upper"]
  ci.lower <- output[ind, "ci.lower"]
  fitmeasures <- fitMeasures(WMH_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi <- fitmeasures["cfi"]
  rmsea <- fitmeasures["rmsea"]
  srmr <- fitmeasures["srmr"]
  tli <- fitmeasures["tli"]
  pheno <- output[ind, "lhs"]
  
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
  results[i,12] <- pheno
  
  Beta_slope= output[snd, "est.std"] 
  SE_slope <- output[snd, "se"]
  p_slope <- output[snd, "pvalue"]
  n_slope <- nobs(WMH_fit_reg)
  ci.upper_slope <- output[snd, "ci.upper"]
  ci.lower_slope <- output[snd, "ci.lower"]
  fitmeasures_slope <- fitMeasures(WMH_fit_reg, c("cfi","rmsea","srmr", "tli"))
  cfi_slope <- fitmeasures_slope["cfi"]
  rmsea_slope <- fitmeasures_slope["rmsea"]
  srmr_slope <- fitmeasures_slope["srmr"]
  tli_slope <- fitmeasures_slope["tli"]
  pheno_slope <- output[snd, "lhs"]
  
  results_slope[i,1] <- i
  results_slope[i,2] <- n_slope
  results_slope[i,3] <- Beta_slope
  results_slope[i,4] <- SE_slope
  results_slope[i,5] <- p_slope
  results_slope[i,6] <- ci.upper_slope
  results_slope[i,7] <- ci.lower_slope
  results_slope[i,8] <- cfi_slope
  results_slope[i,9] <- rmsea_slope
  results_slope[i,10] <- srmr_slope
  results_slope[i,11] <- tli_slope
  results_slope[i, 12] <- pheno_slope
  
}


write.csv(results, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/WMH_intercept_LBC1937_assocs_2_EpiScores_mid_16122022.csv")
write.csv(results_slope, file = "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/WMH_slope_LBC1937_assocs_2_EpiScores_mid_16122022.csv")

