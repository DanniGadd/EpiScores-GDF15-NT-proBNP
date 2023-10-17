#######################################
#### Epigenomic prediction of CRP #####
#######################################

#############################################################
######### Prepare European genetic data for CRP #############
#############################################################

## p17 

## Set working directory 
setwd("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/")

## Load requisite libraries 
library(data.table)


####################################################################################################
#### The format of SNPs (e.g. rsid, CHR:SNP:POS, CHR_SNP_POS) differs a lot between GWAS     #######
#### and it will also proably differ with a target genotype file in a given cohort           #######
#### In order for any PRS software to run, it will need to be able to link up SNPs between   #######
#### your base file (GWAS) and your target file (cohort genotype data). The next step is     #######
#### designed to aid with that. You will need to be sure that the builds are the same too    #######
#### (e.g. hg37 or hg38), as the different SNPs may tag different locations across builds    #######
####################################################################################################

### ROB EXAMPLE FORMATTING

## Read in GWAS data - I downloaded from GWAS Catalog 
# gen1=as.data.frame(fread("/home/robert/GCST90029070_buildGRCh37.tsv"))

####### Generation Scotland formatting #######

## Create column that suits format for Generation Scotland - this is HRC imputed data and has a characteristic SNP format 
gen1$snp=paste(paste(paste(gen1$chromosome, gen1$base_pair_location, sep= "_"),gen1$effect_allele,sep="_"),gen1$other_allele, sep="_")
## Of note, a number of effect alleles (A1) may be the 'other allele' in your genotype file, you can take extra steps to align these 
## but because I worked with multiple cohorts, it wasn't so practical. If working with one cohort (see PREVENT example), you can take extra steps to align base and target file
## Convert p values to numeric format 
gen1$p_value=as.numeric(gen1$p_value)
## Subset to smaller number of SNPs for convenience of analyses
# Here you can normally include all SNPs, commonly all SNPs (P<1) are used to build the PRS but I am using this threshold for consistency with previous CRP studies
gen1=gen1[which(gen1$p_value <= 5e-8), ]
## Save out file for Generation Scotland analyses 
gen1$variant_id=NULL 
#fwrite(gen1, "/Cluster_Filespace/Marioni_Group/Rob/CRP/PGRS/crp_gs_snps.txt", row.names = F, quote = F, sep = "\t")

####### LBC formatting #######

## Create column that suits 1000G format for LBC
gen1$snp=paste(paste(paste(gen1$chromosome, gen1$base_pair_location, sep= ":"),gen1$effect_allele,sep=":"),gen1$other_allele, sep=":")
## Save out file for Generation Scotland analyses 
#fwrite(gen1, "/Cluster_Filespace/Marioni_Group/Rob/CRP/PGRS/crp_lbc_snps.txt", row.names = F, quote = F, sep = "\t")

### You can skip these and these may be found to be sub-optimal. I couldn't see a version of the imputed LBC data in binary files (bed, bim, fam) when I ran this script initially
### I did find the imputed data in VCF format so I changed to binary files and took some extra QC steps 

########################################################################################

### MY ANALYSES FORMATTING

# Load in the GWAS sum stats available for each protein 
GDF <- fread("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GS20K_gdf15_rnk_My8MHH5F8V_imp.stats.gz")
BNP <- fread("/Cluster_Filespace/Marioni_Group/Danni/GS_added_biomarkers/GWAS_caroline/GS20K_nt.probnp_rnk_1ePbnCoIl6_imp.stats.gz")

# GDF15 GS

GDF15 <- GDF[,c('CHR', 'BP', 'ALLELE1', 'ALLELE0', 'BETA', 'SE', 'P_BOLT_LMM', 'SNP')]
names(GDF15) <- c('chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'standard_error', 'p_value', 'snp')
fwrite(GDF15, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gdf_gs_snps.txt", row.names = F, quote = F, sep = "\t")

# GDF15 LBC

GDF15 <- GDF[,c('CHR', 'BP', 'ALLELE1', 'ALLELE0', 'BETA', 'SE', 'P_BOLT_LMM', 'SNP')]
names(GDF15) <- c('chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'standard_error', 'p_value', 'snp')
GDF15$snp <- gsub('_', ':', GDF15$snp)
fwrite(GDF15, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/gdf_lbc_snps.txt", row.names = F, quote = F, sep = "\t")

# BNP GS

BNP1 <- BNP[,c('CHR', 'BP', 'ALLELE1', 'ALLELE0', 'BETA', 'SE', 'P_BOLT_LMM', 'SNP')]
names(BNP1) <- c('chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'standard_error', 'p_value', 'snp')
fwrite(BNP1, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/bnp_gs_snps.txt", row.names = F, quote = F, sep = "\t")

# BNP LBC

BNP2 <- BNP[,c('CHR', 'BP', 'ALLELE1', 'ALLELE0', 'BETA', 'SE', 'P_BOLT_LMM', 'SNP')]
names(BNP2) <- c('chromosome', 'base_pair_location', 'effect_allele', 'other_allele', 'beta', 'standard_error', 'p_value', 'snp')
BNP2$snp <- gsub('_', ':', BNP2$snp)
fwrite(BNP2, "/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/bnp_lbc_snps.txt", row.names = F, quote = F, sep = "\t")


########################################################################################

##### CONVERT LBC 1000G VCF DATA TO 1000G binary files #############
## LBC1936 first 
screen -S change
cd /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/LBC1936_Data/
for x in {1..22} # looping through each chromosome 
do
plink2 --vcf /GWAS_Source/LBC1936/1000G_Phase3v5/LBC1936_1000G_phase3v5_VCF/chr${x}.dose.vcf.gz --make-bed --out /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/LBC1936_Data/chr${x}
--make-bed \
--out chr${x} 
# Here I am finding duplicate SNPs and kicking them out - this was ignored in the LBC data I found so I wanted to make it cleaner 
plink19 --bfile /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/LBC1936_Data/chr${x} --list-duplicate-vars suppress-first --out /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/LBC1936_Data/chr${x} 
plink19 --bfile /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/LBC1936_Data/chr${x} --exclude /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/LBC1936_Data/chr${x}.dupvar --make-bed --out /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_PRS/LBC1936_Data/chr${x}_fix
done

########################################################################################

########### CALCUALTE POLYGENIC RISK SCORES ##########

#### Generation Scotland GDF15 and NT-proBNP scores #########

## You will need to tell it what column names to look for in your formatted base file (GWAS file)
## Binary-target being False means that the phenotype is not binary
## The clumping parameters are taken from previous studies in the group 
## no-regress means that I am not including a phenotype file here. You can get PRSice to do protein ~ PRS for you in your test set to get the R2
## stat beta means that I am using beta and standard errors rather than something else like Z scores or Odds ratios
## I have set a random seed for reproducibility 
## type bed means that I have bed, bim, fam files - other options are available 
## I have only told it to build PRS using SNPs < 5e-8, you can simply add more in here like 1, 0.1, 0.05, 0.01 and so on. 
## Normally people see what one works best in their test set and take that forward, which can be either viewed as pragmatic or cherry-picking. P<5e-8 and P<1 are the most standard/defensible - either using all or just genome-wide significant. I used genome-wide significant due to protein data and wanting to capture 'good' pQTLs
## Of note, plink can run PRS so it can be worth having a play with their documentation. I don't have that to hand but it is often more straightforward than PRSice, other softwares are used too but PRSice is our go to 

## Note that # symbol for the chromosome number, this will essentially run it as a loop and calculate all chromosomes together 
## You will have different output files but an all score file is the PRS at the different thresholds you elected. There will be a log file too reporting the 'bad' SNPs kicked out of analyses

cd /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/

PRSice_linux \
    --A1 effect_allele \
    --A2 other_allele \
    --all \
    --bar-levels 5e-8 \
    --base 00_PRS/gdf_gs_snps.txt \
    --binary-target F \
    --bp base_pair_location \
    --beta \
    --chr chromosome \
    --clump-kb 250 \
    --clump-p 1.000000 \
    --clump-r2 0.25 \
    --fastscore  \
    --model add \
    --no-regress \
    --out 00_PRS/gs_genomewide_gdf \
    --pvalue p_value \
    --seed 817590212 \
    --stat beta \
    --snp snp \
    --se standard_error \
    --target /Cluster_Filespace/Marioni_Group/GS/GS_GWAS/HRC_imputed/GS20K_chr#_HRC.r1-1_nomono_I4_cpra \
    --thread 1 \
    --type bed


cd /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/

PRSice_linux \
    --A1 effect_allele \
    --A2 other_allele \
    --all \
    --bar-levels 5e-8 \
    --base 00_PRS/bnp_gs_snps.txt \
    --binary-target F \
    --bp base_pair_location \
    --beta \
    --chr chromosome \
    --clump-kb 250 \
    --clump-p 1.000000 \
    --clump-r2 0.25 \
    --fastscore  \
    --model add \
    --no-regress \
    --out 00_PRS/gs_genomewide_bnp \
    --pvalue p_value \
    --seed 817590212 \
    --stat beta \
    --snp snp \
    --se standard_error \
    --target /Cluster_Filespace/Marioni_Group/GS/GS_GWAS/HRC_imputed/GS20K_chr#_HRC.r1-1_nomono_I4_cpra \
    --thread 1 \
    --type bed



#### LBC1936 score #########

cd /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/

PRSice_linux \
    --A1 effect_allele \
    --A2 other_allele \
    --all \
    --bar-levels 5e-8 \
    --base 00_PRS/gdf_lbc_snps.txt \
    --binary-target F \
    --bp base_pair_location \
    --beta \
    --chr chromosome \
    --clump-kb 250 \
    --clump-p 1.000000 \
    --clump-r2 0.25 \
    --fastscore  \
    --model add \
    --no-regress \
    --out 00_PRS/lbc1936_genomewide_gdf \
    --pvalue p_value \
    --seed 817590212 \
    --stat beta \
    --snp snp \
    --se standard_error \
    --target 00_PRS/LBC1936_Data/chr#_fix \
    --thread 1 \
    --type bed


#### LBC1936 score #########

cd /Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/

PRSice_linux \
    --A1 effect_allele \
    --A2 other_allele \
    --all \
    --bar-levels 5e-8 \
    --base 00_PRS/bnp_lbc_snps.txt \
    --binary-target F \
    --bp base_pair_location \
    --beta \
    --chr chromosome \
    --clump-kb 250 \
    --clump-p 1.000000 \
    --clump-r2 0.25 \
    --fastscore  \
    --model add \
    --no-regress \
    --out 00_PRS/lbc1936_genomewide_bnp \
    --pvalue p_value \
    --seed 817590212 \
    --stat beta \
    --snp snp \
    --se standard_error \
    --target 00_PRS/LBC1936_Data/chr#_fix \
    --thread 1 \
    --type bed



########################################################################################

