###################################################################################

### LBC file integration 

###################################################################################

# Mean ages at each wave 

# Mean age W5 with SD
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "5") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 

mean(dat3$age, na.rm = T)
sd(dat3$age, na.rm = T)

# Mean age W4 with SD
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "4") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 

mean(dat3$age, na.rm = T)
sd(dat3$age, na.rm = T)

# Mean age W3 with SD
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "3") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 

mean(dat3$age, na.rm = T)
sd(dat3$age, na.rm = T)

# Mean age W2 with SD
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "2") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 

mean(dat3$age, na.rm = T)
sd(dat3$age, na.rm = T)

# Mean age W1 with SD
target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target %>% filter(WAVE == "1") # 801 people
dat3  <- dat3 %>% filter(cohort == "LBC36") # 801 people 

mean(dat3$age, na.rm = T)
sd(dat3$age, na.rm = T)

###################################################################################

### Test EpiScores in W4 

###################################################################################

cd /Local_Scratch/Danni/GDFBNP/04_LBC/20k_GS_training/

screen

R

library(tidyverse)

W4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/00_Updated_results/Results_scores/projections_LBC_W4_241122.csv")

BNP_w4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/LBC_testing/BNP_w4.csv")
GDF_w4 <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/LBC_testing/GDF_w4.csv")


### GDF W4
names(W4)[3] <- 'Basename'
join <- left_join(GDF_w4, W4, by = 'Basename')

library(bestNormalize)
for(i in colnames(join)[c(18)]){ 
  join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
}

cor.test(join$GDF15, join$GDF15_W4) # 0.36

join$age <- join$agedays_w4.x / 365.25
join$sex <- join$sex.x

null <- summary(lm(GDF15_W4 ~ age + sex, data=join))$r.squared
full <- summary(lm(GDF15_W4 ~ age + sex + join$GDF15, data=join))$r.squared
print(round(100*(full - null), 3))

# GDF15 450k
# LBC W4: 8.89%

# Add in PRS for GDF in LBC
prs <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/lbc1936_genomewide_gdf.all.score", header = T)

names(prs)[2] <- 'lbc36no'
prs <- prs[-1]
names(prs)[2] <- 'GDFPRS'

join <- left_join(join, prs, by = 'lbc36no')

null <- summary(lm(GDF15_W4 ~ age + sex, data=join))$r.squared
null_prs <- summary(lm(GDF15_W4 ~ age + sex + GDFPRS, data=join))$r.squared
full <- summary(lm(GDF15_W4 ~ age + sex + join$GDF15, data=join))$r.squared
all <- summary(lm(GDF15_W4 ~ age + sex + join$GDF15 + GDFPRS, data=join))$r.squared
print(round(100*(full - null), 3))

# > null
# [1] 0.09482718
# > null_prs
# [1] 0.1240604
# > full
# [1] 0.1837803
# > all
# [1] 0.23148


table(is.na(join$GDF15))
table(is.na(join$GDF15_W4)) # 322 with both available 

library(ggpubr)
plot1 <- ggplot(join, aes(x=GDF15_W4, y=GDF15)) +
geom_point(colour = "tan1", size = 1) +
geom_smooth(method='lm', colour = "tan1") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) +
xlab("GDF15 protein") + ylab("GDF15 EpiScore") + 
theme_classic() # , label.x = 4.5, label.y = 4.4

# prot_model <- lm(Ig ~ age + sex + GDF15_W4, data = join)
# score_model <- lm(Ig ~ age + sex + GDF15.20k.with.450k.array, data = join)

# join <- join[which(join$enl_peri_space_hl_w2 %in% c('1','0')),]
# prot_model <- glm(enl_peri_space_hl_w2 ~ age + sex + GDF15_W4 + ICV_mm3_wX, 
#   family=binomial(link='logit'), data = join)
# score_model <- glm(enl_peri_space_hl_w2 ~ age + sex + GDF15.20k.with.450k.array + ICV_mm3_wX,
# family=binomial(link='logit'), data = join)


# Compare to GDF15 from grimage

LBC_grim <- read.csv('/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/LBC_clock_output_3489.csv')

target <- read.csv("/Cluster_Filespace/Marioni_Group/LBC/LBC_methylation/archive/target_QC_age_sex_date.csv")
dat3 <- target[which(target$WAVE %in% '4'),]
dat3  <- dat3 %>% filter(cohort == "LBC36") 
dat3 <- dat3 %>% select('Basename', 'set', 'ID')

which(LBC_grim$Basename %in% dat3$Basename)

join <- left_join(join, LBC_grim, by = 'Basename')

# null <- summary(lm(GDF15_W4 ~ age + sex, data=join))$r.squared
# full <- summary(lm(GDF15_W4 ~ age + sex + join$GDF15, data=join))$r.squared
# print(round(100*(full - null), 3))

null <- summary(lm(GDF15_W4 ~ age + sex, data=join))$r.squared
full <- summary(lm(GDF15_W4 ~ age + sex + join$DNAmGDF15, data=join))$r.squared
print(round(100*(full - null), 3))


### BNP w4
names(W4)[3] <- 'Basename'
join <- left_join(BNP_w4, W4, by = 'Basename')

join$bld_NT_ProBNP_w4 <- as.numeric(join$bld_NT_ProBNP_w4)

library(bestNormalize)
for(i in colnames(join)[c(20)]){ 
  join[,i]<- orderNorm(join[,i])$x.t # Rank-Inverse Based Normaliation
}

cor.test(join$NT.proBNP, join$bld_NT_ProBNP_w4) # 0.25

join$age <- join$agedays_w4.x / 365.25
join$sex <- join$sex.x

null <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex, data=join))$r.squared
full <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex + join$NT.proBNP, data=join))$r.squared
print(round(100*(full - null), 3))

# 8.13

# Add in PRS for BNP in LBC
prs <- read.table("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/lbc1936_genomewide_bnp.all.score", header = T)

names(prs)[2] <- 'lbc36no'
prs <- prs[-1]
names(prs)[2] <- 'BNPPRS'

join <- left_join(join, prs, by = 'lbc36no')

null <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex, data=join))$r.squared
null_prs <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex + BNPPRS, data=join))$r.squared
full <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex + join$NT.proBNP, data=join))$r.squared
all <- summary(lm(bld_NT_ProBNP_w4 ~ age + sex + join$NT.proBNP + BNPPRS, data=join))$r.squared
print(round(100*(full - null), 3))

# > null
# [1] 0.007537335
# > null_prs
# [1] 0.04060833
# > full
# [1] 0.08884607
# > all
# [1] 0.09824992

table(is.na(join$NT.proBNP))
table(is.na(join$bld_NT_ProBNP_w4)) #  500 (2 missing values)


library(ggpubr)
plot2 <- ggplot(join, aes(x= bld_NT_ProBNP_w4, y= NT.proBNP)) +
geom_point(colour = "firebrick1", size = 1) +
geom_smooth(method='lm', colour = "firebrick1") + 
theme(axis.text.x=element_text(size=20),     
      axis.text.y=element_text(size=20),
      axis.title.x=element_text(size=20),
      axis.title.y=element_text(size=20)) + 
stat_cor(aes(label = ..r.label..), method = "pearson", cor.coef.name = "r", size = 7.5) +
xlab("Nt-proBNP protein") + ylab("Nt-proBNP EpiScore") + 
theme_classic() # , label.x = 4.5, label.y = 4.4


### PLOT TOGETHER
library(patchwork)
pdf('/Local_Scratch/Danni/GDFBNP/00_Updated_results/Results_scores/CORR_PLOT_JOINT_LBC.pdf', width = 12, height = 5)
plot1 + plot2 
dev.off()

