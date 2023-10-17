
## Extract results for BNP

cd /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/BNP/Outputs/
for i in *_resid.csv

do
 
sigma1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "sigma" | cut -f 1 |  sed 's/:/\n/g' | awk 'NR==1') 
sigma2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "sigma" | cut -f 1 |  sed 's/:/\n/g' | awk 'END{print $NF}') 

A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)

cat $i | cut -d ',' -f $sigma1-$sigma2 > ../Sigma/${B}_output_BNP.csv


beta1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "beta" | cut -f 1 | sed 's/:/\n/g' | awk 'NR==1')
beta2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "beta" | cut -f 1 | sed 's/:/\n/g' | awk 'END{print $NF}')

A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)

cat $i | cut -d ',' -f $beta1-$beta2 > ../Beta/${B}_output_BNP.csv


comp1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "comp" | cut -f 1 | sed 's/:/\n/g' | awk 'NR==1')
comp2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "comp" | cut -f 1 | sed 's/:/\n/g' | awk 'END{print $NF}')

A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)

cat $i | cut -d ',' -f $comp1-$comp2 > ../Comp/${B}_output_BNP.csv

done 

## Extract results for GDF

cd /Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Outputs/
for i in *_resid.csv

do
 
sigma1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "sigma" | cut -f 1 |  sed 's/:/\n/g' | awk 'NR==1') 
sigma2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "sigma" | cut -f 1 |  sed 's/:/\n/g' | awk 'END{print $NF}') 

A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)

cat $i | cut -d ',' -f $sigma1-$sigma2 > ../Sigma/${B}_output_GDF.csv


beta1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "beta" | cut -f 1 | sed 's/:/\n/g' | awk 'NR==1')
beta2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "beta" | cut -f 1 | sed 's/:/\n/g' | awk 'END{print $NF}')

A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)

cat $i | cut -d ',' -f $beta1-$beta2 > ../Beta/${B}_output_GDF.csv


comp1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "comp" | cut -f 1 | sed 's/:/\n/g' | awk 'NR==1')
comp2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "comp" | cut -f 1 | sed 's/:/\n/g' | awk 'END{print $NF}')

A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)

cat $i | cut -d ',' -f $comp1-$comp2 > ../Comp/${B}_output_GDF.csv

done 






####################################################################################################

### PROCESSING THE OUTPUTS FROM BAYESR+

####################################################################################################

### Code taken from daniels cog EWAS repo on gitlab and adapted 

### BNP

## Open R 
R
setwd("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/BNP/Outputs/") 
library(data.table) 

loop = list.files(, pattern = "_resid.csv") 
## Step 1 - Calculate Mean Variance explained by all probes and credible intervals

names <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/BNP/bnp_cpg_subset_list.csv")
names = names$Marker

# for(i in loop){  
i <- loop[1]

  output = matrix(nrow =1, ncol = 1) 
  output <- as.data.frame(output) 
  names(output)[1] <- "Biomarker" 
  
  sigma <- read.csv(paste("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/BNP/Sigma/bnp_output_BNP.csv"))  
  output$Mean_Variance_Explained = mean(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T)
  output$Variance_LCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.025)
  output$Variance_HCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.975)
  output$Biomarker <- "Nt-proBNP"
  write.csv(output, paste0("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/BNP/Processed_results/BNP_sigma_mean_variance_output.csv"), row.names = F) 


### Posterior Inclusion Probability of CpGs 

# Calculate Posterior Inclusion Probability of CpGs
# Over 95% inclusion means significant i.e.  PIP > 0.95 

# Set working directory and file name if looping 

loop = list.files("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/BNP/Comp/", pattern = ".csv") 

# Get the comp file read in as data frame 
  comp <- fread(paste("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/BNP/Comp/bnp_output_BNP.csv"))  
  comp<-as.data.frame(comp) 
dim(comp) # 1000 772667 - so this is 1000 runs kept, for 772667 cpg sites 

# Calculate PIP for each cpg - returns as a list of length 772667
pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})

# Unlist the PIPs so that you have a dataframe of dim  459309 and 1 column with pip for each cpg site 
pip <- as.data.frame(reshape2::melt(unlist(pip)))


names <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/BNP/bnp_cpg_subset_list.csv")
names = names$Marker

# setDT converts to table format 
pip <- setDT(pip, keep.rownames = TRUE) 

# Assign new colnames
names(pip) <- c("Marker", "PIP") 

# Input the order of cpgs saved out in the EWAS step as the markers 
pip$Marker <- names

# This can be looped if more than one biomarker of interest 
pip$Biomarker <- 'Nt-proBNP'

# Order by PIP 
pip2 <- pip[order(-pip$PIP),]

# Which CpGs had PIP greater than 0 
pip3 <- pip2[which(pip2$PIP > 0),] # 421831


# get top hits 
pip4 <- pip2[which(pip2$PIP >= 0.95),]

# Write out
write.csv(pip4, file = paste("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/BNP_PIP_top.csv")) 




### GDF 

## Open R 
R
setwd("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Outputs/") 
library(data.table) 

loop = list.files(, pattern = "_resid.csv") 
## Step 1 - Calculate Mean Variance explained by all probes and credible intervals

names <- read.csv("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/CpGs_GDF/gdf_cpg_subset_list.csv")
names = names$Marker

# for(i in loop){  
i <- loop[1]

  output = matrix(nrow =1, ncol = 1) 
  output <- as.data.frame(output) 
  names(output)[1] <- "Biomarker" 
  
  sigma <- read.csv(paste("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Sigma/gdf_output_GDF.csv"))  
  output$Mean_Variance_Explained = mean(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T)
  output$Variance_LCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.025)
  output$Variance_HCI = quantile(sigma[,2]/rowSums(sigma[,1:2]),na.rm = T, prob =0.975)
  output$Biomarker <- "GDF15"
  write.csv(output, paste0("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/BNP/Processed_results/BNP_sigma_mean_variance_output.csv"), row.names = F) 


### Posterior Inclusion Probability of CpGs 

# Get the comp file read in as data frame 
  comp <- fread(paste("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Comp/gdf_output_GDF.csv"))  
  comp<-as.data.frame(comp) 
dim(comp) # 1000 752722 - so this is 1000 runs kept, for each cpg site

# Calculate PIP for each cpg - returns as a list of length 772667
pip <- lapply(comp,function(x){length(x[which(x>0)])/length(x)})

# Unlist the PIPs so that you have a dataframe of dim  459309 and 1 column with pip for each cpg site 
pip <- as.data.frame(reshape2::melt(unlist(pip)))


names <- read.csv("/Local_Scratch/Danni/GDFBNP/05_MWAS/Preps/CpGs_GDF/gdf_cpg_subset_list.csv")
names = names$Marker

# setDT converts to table format 
pip <- setDT(pip, keep.rownames = TRUE) 

# Assign new colnames
names(pip) <- c("Marker", "PIP") 

# Input the order of cpgs saved out in the EWAS step as the markers 
pip$Marker <- names

# This can be looped if more than one biomarker of interest 
pip$Biomarker <- 'GDF'

# Order by PIP 
pip2 <- pip[order(-pip$PIP),]


# get top hits 
pip4 <- pip2[which(pip2$PIP >= 0.95),]

# Write out
write.csv(pip4, file = paste("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/GDF_PIP_top.csv")) 


####################################################################################################

### Get the summary table for the CpGs for table 1

gdf <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/GDF_PIP_top.csv")
bnp <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/BNP_PIP_top.csv")

t <- rbind(gdf, bnp)
write.csv(t, '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/joint.csv', row.names = F)


res <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/results.txt")
stud <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/studies.txt")
res <- merge(res, stud, by = 'StudyID')

EWAS <- res
LBCj <- t

# Now look at how many proteins for each of the cpgs and consider annotating them to EWAS catalog 

names(LBCj)[2] <- "CpG_site"
library(tidyverse)
counts3 <- count(LBCj, CpG_site)

results <- counts3

# Set a place for the annotations 
results$Gene <- "X"
results$Annotation <- "X"

list <- list()

# Get annotations for the 14 cpgs of interest 

for (i in 1:14){
	print(i)
	cpg <- results[i,1]
	anno <- EWAS
	anno_cpg <- anno[which(anno$CpG %in% cpg),]
	anno_cpg <- anno_cpg[which(anno_cpg$P < 3.6e-8),]
	anno_cpg <- anno_cpg[which(anno_cpg$N >= 100),]
	trait <- anno_cpg$Trait %>% unique()
	str <- str_c(trait, collapse = ", ")
	results[i,3] <- unique(anno_cpg$Gene)
	results[i,4] <- str

	list[[i]] <- anno_cpg
	
}

list2 <- do.call(rbind, list)

# Save off results file 
write.csv(results, "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/suppl_table_cpgs_counted_and_annotated.csv", row.names = F)

# Save off full set of indexing for all cpgs joint for suppl
write.csv(list2, "/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/suppl_table_all_cpgs_info.csv", row.names = F)


##################################################################################

# Check for previous EWAS of each protein

res <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/results.txt")
stud <- read.delim("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/studies.txt")
res <- merge(res, stud, by = 'StudyID')

test <- res[grep('probnp', res$Trait),]
test <- res[grep('proBNP', res$Trait),]
test <- res[grep('BNP', res$Trait),]

test <- res[grep('GDF15', res$Trait),]

write.csv(test, '/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/GDF_previous.csv', row.names = F)

cpgs <- read.csv('/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/00_MWAS_20k/Results/GDF_PIP_top.csv')

which(cpgs$Marker %in% test$CpG)