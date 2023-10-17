####################################################################################

### ELNETS for the 2 new biomarkers in GS to train DNAm episcores 

####################################################################################

cd /Local_Scratch/Danni/GDFBNP/

screen 

R

# Load packages and functions for use in models 
library(bigmemory)
library(optparse)
library(biglasso)
library(tidyverse)
library(foreign)

# EPIC files first, then will train 450k option
BNP_train <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/ALL_BNP.csv")
GDF_train <- read.csv("/Cluster_Filespace/Marioni_Group/Danni/EpiScores_20k/GDFBNP/00_All_initial_files/02_Relatedness_mapping/ALL_GDF15.csv")


################################################################################

# Split big lasso matricies using cbindBM - load relevant functions for this

################################################################################

# x needs to be a list of big matrices, and they have to have the same rownames
cbindBM_list <- function(x, binding="right", 
                         z=NULL, type=NULL, separated=NULL,
                         backingfile=NULL, backingpath=NULL,
                         descriptorfile=NULL, binarydescriptor=FALSE,
                         shared=TRUE, erase = TRUE)
{
  
  if (is.null(type)) type <- typeof(x[[1]])
  if (is.big.matrix(x[[1]])) {
    if (is.null(separated)) separated <- is.separated(x[[1]])
  } else {
    separated <- FALSE
  }
  
  cols_list <- list()
  total_cols <- 0
  for (i in 1:length(x)) {
    cols <- cleanupcols(NULL, ncol(x[[i]]), colnames(x[[i]]))
    cols_list <- append(cols_list, list(cols))
    total_cols <- total_cols + ncol(x[[i]])
  }    
  
  if (is.null(z)) {
    z <- big.matrix(nrow=nrow(x[[1]]), ncol=total_cols, type=type, init=NULL,
                    dimnames=dimnames(x[[1]]), separated=separated,
                    backingfile=backingfile, backingpath=backingpath,
                    descriptorfile=descriptorfile,
                    binarydescriptor=binarydescriptor, shared=shared)
  }
  
  counter <- 0
  for (i in 1:length(cols_list)) {
    print(i)
    if (i == 1) {
      z[, 1:length(cols_list[[i]])] <- x[[i]][,cols_list[[i]]]
    } else {
      z[, (counter + 1):(counter + length(cols_list[[i]]))] <- x[[i]][,cols_list[[i]]]
    }
    counter <- counter + length(cols_list[[i]])
    print(counter)
    
    if (erase == TRUE) {
      cat("\nErasing chunk and liberating memory...\n\n")
      x[[i]] <- "Replacement"
      gc()
    }
  }
  return(z)
}


cleanupcols <- function(cols=NULL, nc=NULL, colnames=NULL) {
  if (is.null(cols)) cols <- 1:nc
  else {
    if (!is.numeric(cols) & !is.character(cols) & !is.logical(cols))
      stop("column indices must be numeric, logical, or character vectors.")
    if (is.character(cols))
      if (is.null(colnames)) stop("column names do not exist.")
    else cols <- mmap(cols, colnames)
    if (is.logical(cols)) {
      if (length(cols) != nc)
        stop(paste("column vector length must match the number of",
                   "columns of the matrix."))
      cols <- which(cols)
    }
    tempj <- .Call("CCleanIndices", as.double(cols), as.double(nc), PACKAGE="bigmemory")
    if (is.null(tempj[[1]])) stop("Illegal column index usage in extraction.\n")
    if (tempj[[1]]) cols <- tempj[[2]]
  }
  return(cols)
}


####################################################################################

### 20k, GDF, 450k, 20cv, pQTL unadjusted

####################################################################################

# 20k x file
all <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/20k_450K.rds")

GDF_x_train <- all[which(rownames(all) %in% GDF_train$Sample_Sentrix_ID),]

GDF_x_train <- GDF_x_train[match(GDF_train$Sample_Sentrix_ID, rownames(GDF_x_train)),]
identical(GDF_train$Sample_Sentrix_ID, rownames(GDF_x_train))

meth <- GDF_x_train
div <- 20 # Number of chunks to divide OG methylation dataframe

por <- ceiling(length(colnames(meth))/div)
chunk_list <- list()

for (i in 1:div) {
  cat(paste0("\nWorking on chunk: ", i, " of ", div))
  if (i == 1) {
    chunk <- as.big.matrix(meth[,1:(por-1)])
  } else if (i == div) {
    chunk <- as.big.matrix(meth[,(por*(i-1)):length(colnames(meth))])
  } else {
    chunk <- as.big.matrix(meth[,(por*(i-1)):((por*i)-1)])
  }
  cat("\nMade chunk. Appending to chunk list...\n")
  chunk_list <- append(chunk_list, list(chunk))
  gc()
}

# Saving names prior to chunk fusing
names <- colnames(meth)
rm(meth)

cat("\nRAM clean up...\n\n")
gc()

cat("\nFusing chunks!\n\n")
GDF_x <- cbindBM_list(x = chunk_list)
rm(chunk, chunk_list)

# Set CpG names
options(bigmemory.allow.dimnames=TRUE)
colnames(GDF_x)<- names



################################################################################

### RUN

################################################################################

location <- "/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/"

set.seed(1783) # set seed to ensure fold variation minimised 

# Set list of traits to run through 
list <- names(GDF_train) # we want to run 2 and 3 as columns with/without pQTL adjustments

# Check identical order of individuals
identical(rownames(GDF_x), GDF_train$Sample_Sentrix_ID) # TRUE

# Assign ytrain
ytrain <- GDF_train

i <- 2
  q <- ytrain[1] # Get just basenames for people in the y variable 
  p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable 
  name_p <- colnames(p) # Get the name of the protein for this iteration
  y <- cbind(q,p) # Bind Basename and protein data together into one set 
  names(y)[2] <- "pheno" # Assign a generic name to the protein variable
  y <- as.numeric(y$pheno) # Create a numeric list for this variable to feed in as y to the model

  # Run training 
  lasso.cv <- cv.biglasso(GDF_x, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
  fit <- biglasso(GDF_x, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit 
  coefs <- coef(fit) # Extract coeficients 
  coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
  coefs$Predictor <- 'GDF15' # Assign protein identifier
  names(coefs)[1] <- "Coefficient" # Tidy naming 
  coefs$CpG <- rownames(coefs) # Create episcores column
  coefs <- coefs[-1,]
  coefs <- coefs[c(2,3,1)]

  # Save coefficients
  write.csv(coefs, file = paste0(location, name_p, "20k_450K_20cv_pQTL_unadjusted.csv"), row.names = F) 

################################################################################

### 20k, BNP, 450k, 20cv

################################################################################

# 20k x file
all <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/20k_450K.rds")

BNP_x_train <- all[which(rownames(all) %in% BNP_train$Sample_Sentrix_ID),]

all = NULL

BNP_x_train <- BNP_x_train[match(BNP_train$Sample_Sentrix_ID, rownames(BNP_x_train)),]
identical(BNP_train$Sample_Sentrix_ID, rownames(BNP_x_train))

meth <- BNP_x_train

div <- 15 # Number of chunks to divide OG methylation dataframe

por <- ceiling(length(colnames(meth))/div)
chunk_list <- list()

for (i in 1:div) {
  cat(paste0("\nWorking on chunk: ", i, " of ", div))
  if (i == 1) {
    chunk <- as.big.matrix(meth[,1:(por-1)])
  } else if (i == div) {
    chunk <- as.big.matrix(meth[,(por*(i-1)):length(colnames(meth))])
  } else {
    chunk <- as.big.matrix(meth[,(por*(i-1)):((por*i)-1)])
  }
  cat("\nMade chunk. Appending to chunk list...\n")
  chunk_list <- append(chunk_list, list(chunk))
  gc()
}

# Saving names prior to chunk fusing
names <- colnames(meth)
rm(meth)

cat("\nRAM clean up...\n\n")
gc()

cat("\nFusing chunks!\n\n")
BNP_x <- cbindBM_list(x = chunk_list)
rm(chunk, chunk_list)

# Set CpG names
options(bigmemory.allow.dimnames=TRUE)
colnames(BNP_x)<- names



################################################################################

### RUN

################################################################################

location <- "/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/"

set.seed(1783) # set seed to ensure fold variation minimised 

# Set list of traits to run through 
list <- names(BNP_train)

# Check identical order of individuals
identical(rownames(BNP_x), BNP_train$Sample_Sentrix_ID) # TRUE

# Assign ytrain
ytrain <- BNP_train

  i <- 2 # set due to y varibales only having 2 cols 
  q <- ytrain[1] # Get just basenames for people in the y variable 
  p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable 
  name_p <- colnames(p) # Get the name of the protein for this iteration
  y <- cbind(q,p) # Bind Basename and protein data together into one set 
  names(y)[2] <- "pheno" # Assign a generic name to the protein variable
  y <- as.numeric(y$pheno) # Create a numeric list for this variable to feed in as y to the model

  # Run training 
  lasso.cv <- cv.biglasso(BNP_x, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
  fit <- biglasso(BNP_x, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit 
  coefs <- coef(fit) # Extract coeficients 
  coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
  coefs$Predictor <- 'NT-proBNP' # Assign protein identifier
  names(coefs)[1] <- "Coefficient" # Tidy naming 
  coefs$CpG <- rownames(coefs) # Create episcores column
  coefs <- coefs[-1,]
  coefs <- coefs[c(2,3,1)]

  # Save coefficients
  write.csv(coefs, file = paste0(location, name_p, "20k_450K_20cv_pQTL_unadjusted.csv"), row.names = F) 



####################################################################################

### GDF15 EPIC 

####################################################################################

# 20k x file
all <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/20k_EPIC.rds")

GDF_x_train <- all[which(rownames(all) %in% GDF_train$Sample_Sentrix_ID),]

all = NULL

GDF_x_train <- GDF_x_train[match(GDF_train$Sample_Sentrix_ID, rownames(GDF_x_train)),]
identical(GDF_train$Sample_Sentrix_ID, rownames(GDF_x_train))

meth <- GDF_x_train
div <- 20 # Number of chunks to divide OG methylation dataframe

por <- ceiling(length(colnames(meth))/div)
chunk_list <- list()

for (i in 1:div) {
  cat(paste0("\nWorking on chunk: ", i, " of ", div))
  if (i == 1) {
    chunk <- as.big.matrix(meth[,1:(por-1)])
  } else if (i == div) {
    chunk <- as.big.matrix(meth[,(por*(i-1)):length(colnames(meth))])
  } else {
    chunk <- as.big.matrix(meth[,(por*(i-1)):((por*i)-1)])
  }
  cat("\nMade chunk. Appending to chunk list...\n")
  chunk_list <- append(chunk_list, list(chunk))
  gc()
}

# Saving names prior to chunk fusing
names <- colnames(meth)
rm(meth)

cat("\nRAM clean up...\n\n")
gc()

cat("\nFusing chunks!\n\n")
GDF_x <- cbindBM_list(x = chunk_list)
rm(chunk, chunk_list)

# Set CpG names
options(bigmemory.allow.dimnames=TRUE)
colnames(GDF_x)<- names



################################################################################

### W1W3 train - RUN GDF15 TRAINING

################################################################################

location <- "/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/"

set.seed(1783) # set seed to ensure fold variation minimised 

# Set list of traits to run through 
list <- names(GDF_train) # we want to run 2 and 3 as columns with/without pQTL adjustments

# Check identical order of individuals
identical(rownames(GDF_x), GDF_train$Sample_Sentrix_ID) # TRUE

# Assign ytrain
ytrain <- GDF_train

i <- 2
  q <- ytrain[1] # Get just basenames for people in the y variable 
  p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable 
  name_p <- colnames(p) # Get the name of the protein for this iteration
  y <- cbind(q,p) # Bind Basename and protein data together into one set 
  names(y)[2] <- "pheno" # Assign a generic name to the protein variable
  y <- as.numeric(y$pheno) # Create a numeric list for this variable to feed in as y to the model

  # Run training 
  lasso.cv <- cv.biglasso(GDF_x, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
  fit <- biglasso(GDF_x, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit 
  coefs <- coef(fit) # Extract coeficients 
  coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
  coefs$Predictor <- 'GDF15' # Assign protein identifier
  names(coefs)[1] <- "Coefficient" # Tidy naming 
  coefs$CpG <- rownames(coefs) # Create episcores column
  coefs <- coefs[-1,]
  coefs <- coefs[c(2,3,1)]

  # Save coefficients
  write.csv(coefs, file = paste0(location, name_p, "20k_EPIC_20cv.csv"), row.names = F) 

################################################################################

### W1W3 train - PREP BNP DNAm

################################################################################

# 20k x file
all <- readRDS("/Local_Scratch/Danni/GDFBNP/00_Preps/20k_EPIC.rds")

BNP_x_train <- all[which(rownames(all) %in% BNP_train$Sample_Sentrix_ID),]

all = NULL

BNP_x_train <- BNP_x_train[match(BNP_train$Sample_Sentrix_ID, rownames(BNP_x_train)),]
identical(BNP_train$Sample_Sentrix_ID, rownames(BNP_x_train))

meth <- BNP_x_train
div <- 20 # Number of chunks to divide OG methylation dataframe

por <- ceiling(length(colnames(meth))/div)
chunk_list <- list()

for (i in 1:div) {
  cat(paste0("\nWorking on chunk: ", i, " of ", div))
  if (i == 1) {
    chunk <- as.big.matrix(meth[,1:(por-1)])
  } else if (i == div) {
    chunk <- as.big.matrix(meth[,(por*(i-1)):length(colnames(meth))])
  } else {
    chunk <- as.big.matrix(meth[,(por*(i-1)):((por*i)-1)])
  }
  cat("\nMade chunk. Appending to chunk list...\n")
  chunk_list <- append(chunk_list, list(chunk))
  gc()
}

# Saving names prior to chunk fusing
names <- colnames(meth)
rm(meth)

cat("\nRAM clean up...\n\n")
gc()

cat("\nFusing chunks!\n\n")
BNP_x <- cbindBM_list(x = chunk_list)
rm(chunk, chunk_list)

# Set CpG names
options(bigmemory.allow.dimnames=TRUE)
colnames(BNP_x)<- names



################################################################################

### W1W3 train - RUN BNP TRAINING

################################################################################

location <- "/Local_Scratch/Danni/GDFBNP/00_Updated_results/GS_EPIC_pQTL/"

set.seed(1783) # set seed to ensure fold variation minimised 

# Set list of traits to run through 
list <- names(BNP_train)

# Check identical order of individuals
identical(rownames(BNP_x), BNP_train$Sample_Sentrix_ID) # TRUE

# Assign ytrain
ytrain <- BNP_train

  i <- 2 # set due to y varibales only having 2 cols 
  q <- ytrain[1] # Get just basenames for people in the y variable 
  p <- ytrain[i] # Get the protein data for the iteration of interest from the y variable 
  name_p <- colnames(p) # Get the name of the protein for this iteration
  y <- cbind(q,p) # Bind Basename and protein data together into one set 
  names(y)[2] <- "pheno" # Assign a generic name to the protein variable
  y <- as.numeric(y$pheno) # Create a numeric list for this variable to feed in as y to the model

  # Run training 
  lasso.cv <- cv.biglasso(BNP_x, y, family="gaussian", alpha = 0.5, ncores = 8, nfolds = 20) # cross validation to get best lambda
  fit <- biglasso(BNP_x, y, family = "gaussian", alpha = 0.5, ncores = 8, lambda = lasso.cv$lambda.min) # model fit 
  coefs <- coef(fit) # Extract coeficients 
  coefs <- as.data.frame(coefs[which(coefs!=0),]) # Remove coeficients that are 0 
  coefs$Predictor <- 'NT-proBNP' # Assign protein identifier
  names(coefs)[1] <- "Coefficient" # Tidy naming 
  coefs$CpG <- rownames(coefs) # Create episcores column
  coefs <- coefs[-1,]
  coefs <- coefs[c(2,3,1)]

  # Save coefficients
  write.csv(coefs, file = paste0(location, name_p, "W1W3_train_weights_EPIC.csv"), row.names = F) 

# 445 for EPIC 


##########################################################################################################

# /dev/shm
# rm -rf /dev/shm/*
# pkill -u danni

# sessionInfo()

