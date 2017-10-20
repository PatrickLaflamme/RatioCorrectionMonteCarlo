#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#check what the commandArgs are requesting
if(length(args) > 0){
  
  # Should we run with MPI?
  run_mpi <- length(grep("mpi", args))
  
  # Where to save the data?
  path <- grep("--path", args) + 1
  
  if( run_mpi ){
    library(Rmpi)
    mpi.setup.rngstream(12345)
    mpi.bcast.Robj2slave(all=TRUE)
  } 
  
  if( length(path) ) {
    path <- args[path]
  } else {
    path <- '~/Desktop/Test/'
  }
} else {
  path <- '~/Desktop/Test/'
}

library(moments)

# Function to generate the index values from the saved data.
genIndex <- function(data, condition){
  ## Generates an index score for each data point. Generates it in one of 3 ways, as listed in the 'type' of condition.
  ##
  ## data      :: always has shape: sampleSize x 2 x n_samples
  ##                                The 2 columns are: random sample dY and slope
  ## condition :: type 1: participant-wise difference/slope - Index method
  ##              type 2: participant-wise difference/overall slope - Zero variance methon
  ##              type 3: taylor-corrected effect scores based on Feiller's theorem - Feiler's method
  
  if(condition == 1){
    indexData <- data[,1,]/data[,2,] #divide each individual effect by each individual slope
  }
  
  else if(condition ==2){
    indexData <- t(t(data[,1,])/apply(data[,2,], 2, mean)) #divide each individual effect by the sample mean slope
  }
  
  else if (condition==3){
    meanSlope <- apply(data[,2,], 2, mean) #calculate the sample mean slope
    meanEffect <- apply(data[,1,], 2, mean) #calculate the sample mean effect
    
    indexData <- t(matrix(meanEffect/meanSlope, ncol=dim(data)[1], nrow=dim(data)[3])) * (1 + t(t(data[,1,])/meanEffect) - t(t(data[,2,])/meanSlope)) #perform the correction
  }
  else if(condition==4){
    meanSlope <- apply(data[,2,], 2, mean) #calculate the sample mean slope
    meanEffect <- apply(data[,1,], 2, mean) #calculate the sample mean effect
    
    indexData <- data[,1,] - t(matrix(meanEffect/meanSlope, ncol=dim(data)[1], nrow=dim(data)[3])) * data[,2,]
  }
  return(indexData)
}


test_between_datasets <- function(datapaths, condition = 1){
  # Function that takes dataset filepaths as input, and returns:
  # - the skew of the difference
  # - the mean of the difference
  # - the type 1 error rate if no true effect is present
  # - the type 2 error rate if a true effect is present.
  #
  # TODO: extract real effect sizes from file names in order to test if we are looking for alpha or beta values.
  #
  # Inputs:
  # datapaths:    A vector of 2 filepaths for the two datapaths.
  # condition:    Metod to use to generate participants' index scores.
  
  infoA <- strsplit(basename(datapaths[1]),"_")[[1]]
  infoB <- strsplit(basename(datapaths[2]),"_")[[1]]
  
  effectA <- as.numeric(infoA[1])
  effectB <- as.numeric(infoB[1])
  
  effectDiff <- effectA - effectB
  
  # Load and prepare dataA
  dataA <- readRDS(datapaths[1])
  indexA <- genIndex(dataA, condition)
  
  # Load and prepare dataA
  dataB <- readRDS(datapaths[2])
  indexB <- genIndex(dataB, condition)
  
  # Create sampling distributions for indexes
  meanDiffs <- colMeans(indexA) - colMeans(indexB)
  
  # Calculate the skew of the sampling distributions.
  skewDiffs <- skewness(meanDiffs)
  
  kurtosisDiffs <- kurtosis(meanDiffs)
  
  sum_vars <- apply(indexA,2,var) + apply(indexB,2,var)
  
  tValues <- meanDiffs/ (sapply(sum_vars,sqrt)/dim(indexA)[2])
  
  pValues <- sapply(tValues, pt, df = dim(indexA)[2])
  
  
  if(effectDiff==0){
  error <- sum(pValues<0.05)/dim(indexA)[2]
  errorType <- "Alpha"
  } else {
  error <- sum(pValues>=0.05)/dim(indexA)[2]
  errorType <- "Beta"
  }
  
  return(c(effectDiff, skewDiffs, kurtosisDiffs, errorType, error))
}

# List files provided in the 
files <- list.files(path, pattern = ".rds", full.names = T, include.dirs = F)

# get a matrix of all unique comparisons between datasets.
all_comparisons <- combn(files, m = 2)

# If we're running with MPI, apply in parallel, otherwise apply in serial
if(run_mpi){
  output <- mpi.parApply(all_comparisons, 2, test_between_datasets)
} else {
  output <- apply(all_comparisons, 2, test_between_datasets)
}

output <- data.frame(t(output))

colnames(output) <- c("Effect Difference", "Skewness", "Kurtosis", "Error Type", "Error")

dir.create(paste(path,"output", sep='/'), showWarnings = F)

saveRDS(output, file=paste0(path, "output/between_subjects_t.rds"))
