##################################################################################

### LBC1936 - SEM results processing 

##################################################################################

# /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/00_Updated_results/Results_LBC_brain/

screen

R

# Load in packages 
library(readxl)
library(tidyverse)

# Basic model collate files 

path <- '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/'

files <- list.files(path)
files <- files[grep('basic', files)]

res <- list()


for (i in 1:length(files)){
    location <- files[i]
    file <- read.csv(paste0(path, location))
    file <- file[-1]
    res[[i]] <- file
    print(file)
}

basic <- do.call(rbind, res)
basic <- basic[order(basic$P),]
basic$FDR <- p.adjust(basic$P, method = 'BH')
length(which(basic$FDR < 0.05))


write.csv(basic, '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/Results_processed/basic.csv', row.names = F)

# Mid model collate files 

path <- '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/'

files <- list.files(path)
files <- files[grep('mid', files)]

res <- list()


for (i in 1:length(files)){
    location <- files[i]
    file <- read.csv(paste0(path, location))
    file <- file[-1]
    res[[i]] <- file
    print(file)
}

full <- do.call(rbind, res)
full <- full[order(full$P),]
# full$FDR <- p.adjust(full$P, method = 'BH')
length(which(full$P < 0.05))


write.csv(full, '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/Results_processed/full.csv', row.names = F)



# Full model collate files 

path <- '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/'

files <- list.files(path)
files <- files[grep('full', files)]

res <- list()


for (i in 1:length(files)){
    location <- files[i]
    file <- read.csv(paste0(path, location))
    file <- file[-1]
    res[[i]] <- file
    print(file)
}

full <- do.call(rbind, res)
full <- full[order(full$P),]
# full$FDR <- p.adjust(full$P, method = 'BH')
length(which(full$P < 0.05))



write.csv(full, '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/Results_processed/WBC.csv', row.names = F)




#####################################################################################

### PROCESS RESULTS

# Basic = age and sex
# full = basic + lifestyle/health covars
# WBC = full + WBC estimates (i.e. the maximal model)

basic <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/Results_processed/basic.csv")
full <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/Results_processed/full.csv")
WBC <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/Results_processed/WBC.csv")

basic$retain <-  paste0(basic$SeqId, basic$phenotype)
full$retain <-  paste0(full$SeqId, full$phenotype)
# WBC$retain <-  paste0(WBC$SeqId, WBC$phenotype)

keep <- basic[which(basic$FDR < 0.05),]

full_keep <- full[which(full$retain %in% keep$retain),]
full_assoc <- full_keep[which(full_keep$P < 0.05),]

full$assoc <- ifelse(full$retain %in% full_assoc$retain, 'red', 'black')

# WBC_keep <- WBC[which(WBC$retain %in% keep$retain),]
# WBC_assoc <- WBC_keep[which(WBC_keep$P < 0.05),]

join <- left_join(basic, full, by = 'retain')
# join <- left_join(join, WBC, by = 'retain')

write.csv(join, '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/Results_processed/joint.csv', row.names = F)

###############################################################

### PLOT RESULTS

library(tidyverse)
library(ggplot2)

bind <- full_assoc

bind$Protein <- c('GDF15 EpiScore', 'GDF15 EpiScore', 'GDF15 EpiScore',  'NT-proBNP EpiScore', 'NT-proBNP EpiScore')
bind$Phenotype <- c('Normal appearing white matter volume',
    'Total brain volume', 'General cognitive ability', 'Normal appearing white matter volume',
    'Total brain volume')

bind$Naming <- paste0(bind$Protein, ' - ', bind$Phenotype)

# Add colour column assignment
bind = bind %>% mutate(Col = case_when(
  bind$beta < 0 ~ "royalblue",
  bind$beta > 0 ~ "tomato2"))

# Set naming to match ordering 
bind$Naming = factor(bind$Naming, levels=unique(bind$Naming[order(bind$beta)]))

# Set colours to match ordering 
bind$Col = factor(bind$Col, levels=unique(bind$Col[order(bind$beta)]))

My_Theme = theme(
  axis.title.x = element_text(size = 20),
  axis.text.x = element_text(size = 20),
  axis.text.y = element_text(size = 20),
  axis.title.y = element_text(size = 24),
  strip.text = element_text(size = 20, face = "bold"),
  legend.text=element_text(size=24),
  legend.title=element_text(size=24, face = "bold"), legend.position = "none")

pdf("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_Analysis_20k/04_LBC/SEM_results/LBC_brain/Results_processed/BRAIN_JOINT_RESULTS_PLOT.pdf", width = 20, height = 4)
ggplot(bind, aes(x = beta, y = Naming, color = bind$Col)) + 
    geom_point(size = 4.5, color = bind$Col) +
        geom_errorbarh(aes(xmax = ci.upper, xmin = ci.lower), size = .9, height = 
                    .4, color = bind$Col) +
        geom_vline(aes(xintercept = 0), size = .25, linetype = "dashed") + 
    coord_trans(x = scales:::exp_trans(10)) + theme_classic() + theme_classic() + 
    My_Theme + 
    theme(panel.grid.minor = element_blank()) +
    ylab("") +
    xlab("Beta") +
    ggtitle("") + theme(plot.title = element_text(hjust = 0.5, size = 27)) + xlim(-0.4, 0)
dev.off()
