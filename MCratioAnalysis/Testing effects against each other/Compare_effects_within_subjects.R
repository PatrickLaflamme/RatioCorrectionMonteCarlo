#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#Read through for code errors like the t-values error.

#check what the commandArgs are requesting
if(length(args) > 0){
  
  # Should we run with MPI?
  run_mpi <- length(grep("mpi", args))
  
  # Where to save the data?
  path <- grep("--path", args) + 1
  
  # Where to save the data?
  condition <- grep("--condition", args) + 1
  
  if( run_mpi ){
    library(Rmpi)
    .Last <- function(){ 
      if (is.loaded("mpi_initialize")){ 
        if (mpi.comm.size(1) > 0){ 
          print("Please use mpi.close.Rslaves() to close slaves.") 
          mpi.close.Rslaves() 
        } 
        print("Please use mpi.quit() to quit R") 
        .Call("mpi_finalize") 
      } 
    }
  } 
  
  if( length(path) ) {
    path <- args[path]
  } else {
    path <- '~/Desktop/Test/'
  }
  if( length(condition) ) {
    condition <- args[condition]
  } else {
    condition <- 1
  }
} else {
  run_mpi <- 0
  path <- '~/Desktop/Test/'
  condition <- 2
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


test_between_datasets <- function(datapath, condition){
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
  
  info <- strsplit(basename(datapath),"_")[[1]]
  
  effectA <- as.numeric(info[1])
  effectB <- as.numeric(info[5])
  
  varEffectA <- as.numeric(info[2])
  varEffectB <- as.numeric(info[6])
  
  slopeA <- as.numeric(info[3])
  slopeB <- as.numeric(info[7])
  
  varSlopeA <- as.numeric(info[4])
  varSlopeB <- as.numeric(info[8])
  
  HiddenEffectscorrelation <- as.numeric(info[9])
  slopescorrelation <- as.numeric(info[10])
  effectSlopecorrelationA <- as.numeric(info[11])
  effectSlopecorrelationB <- as.numeric(strsplit(info[12], ".rds")[[1]][1])
  
  effectCoeffVarA <- varEffectA/effectA
  effectCoeffVarB <- varEffectB/effectB
  
  slopeCoeffVarA <- varSlopeA/slopeA
  slopeCoeffVarB <- varSlopeB/slopeB
  
  effectDiff <- effectA - effectB
  
  obsEffectA <- effectA/slopeA
  obsEffectB <- effectB/slopeB
  
  obsEffectDiff <- obsEffectA - obsEffectB
  
  # Load and prepare data
  data <- readRDS(datapath)
  indexA <- genIndex(data[,c(1,2),], condition)
  indexB <- genIndex(data[,c(3,4),], condition)
  
  # Create sampling distributions for indexes
  meanDiffs <- colMeans(indexA) - colMeans(indexB)
  
  # Calculate the skew of the sampling distributions.
  skewDiffs <- skewness(meanDiffs)
  
  kurtosisDiffs <- kurtosis(meanDiffs)
  
  sum_vars <- apply(indexA,2,var) + apply(indexB,2,var)
  
  tValues <- meanDiffs/ (sapply(sum_vars/dim(indexA)[1],sqrt))
  
  pValues <- sapply(-abs(tValues), pt, df = dim(indexA)[1] - 2)
  
  
  if(effectDiff==0){
    error <- sum(pValues<0.05)/dim(indexA)[2]
    errorType <- 0
  } else {
    error <- sum(pValues>=0.05)/dim(indexA)[2]
    errorType <- 1
  }
  
  return(c(effectDiff, obsEffectDiff, obsEffectA, effectCoeffVarA, slopeCoeffVarA, obsEffectB, effectCoeffVarB, slopeCoeffVarB,  HiddenEffectscorrelation, slopescorrelation, effectSlopecorrelationA, effectSlopecorrelationB, skewDiffs, kurtosisDiffs, errorType, error))
}

# List files provided in the 
files <- list.files(path, pattern = ".rds", full.names = T, include.dirs = F)

# If we're running with MPI, apply in parallel, otherwise apply in serial
if(run_mpi){
  output <- mpi.parApply(files, 2, test_between_datasets, condition = condition)
} else {
  output <- sapply(files, test_between_datasets, condition = condition, simplify = T, USE.NAMES = F)
}

output <- data.frame(t(output), stringsAsFactors = F)

colnames(output) <- c("Effect Difference","Expected Observed Effect Difference", "Observed Effect 1", "Dataset 1 effect coefficient of variance", "Dataset 1 slope coefficient of variance", "Observed Effect 2", "Dataset 2 effect coefficient of variance","Dataset 2 slope coefficient of variance", "Hidden Effect Correlation", "Slopes correlation", "Hidden Effect-Slope correlation dataset 1", "Hidden Effect-Slope correlation dataset 2", "Skewness", "Kurtosis", "Error Type", "Error")
output$`Error Type` <- factor(output$`Error Type`, labels = c("Alpha","Beta"))

dir.create(paste(path,"output", sep='/'), showWarnings = F)

saveRDS(output, file=paste0(path, "/output/within_subjects_t_condition_", condition, ".rds"))

if(run_mpi){
  mpi.close.Rslaves() 
  mpi.quit()
}
