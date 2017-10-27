#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)


# Should we run with MPI?
run_mpi <- 1

# Where to save the data?
path <- "RatioCorrectionMonteCarlo/MCratioAnalysis/Testing\ effects\ against\ each\ other/Generate_within_subject_data_files.R"

# Where to save the data?
n_samples <- 10000

# get the number of cores
cores <- 10

#condition to run
condition <- 2
 
library(Rmpi)
library(mnormt)

# data generation function
Datagen <- function(n_samples, sample_size, means,varcov, seed){
  
  #create a sample_size x n_samples array x groups
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  sample <- array(rmnorm(n_samples*sample_size,means,varcov), c(sample_size, n_samples, length(means)))
  
  #now change it to a more intelligible sample_size x groups x n_samples array
  sample <- aperm(sample, c(1,3,2))
  
  return(sample)
  
}

# create data and save to file function
gen_data_file <- function(varVector, dirpath, n_samples, sample_size=30, seed = 1234){
  
  HiddenEffectA <- varVector[1]
  HiddenEffectB <- varVector[2]
  EffectVarA <- varVector[3]
  EffectVarB <- varVector[4]
  slopeA <- varVector[5]
  slopeB <- varVector[6]
  slopeVarA <- varVector[7]
  slopeVarB <- varVector[8]
  HiddenEffectscorrelation <- varVector[9]
  slopescorrelation <- varVector[10]
  effectSlopecorrelationA <- varVector[11]
  effectSlopecorrelationB <- varVector[12]
  
  varcov <- diag(c(EffectVarA, slopeVarA,EffectVarB,slopeVarB))
  
  varcov[1,3] <- HiddenEffectscorrelation * sqrt(EffectVarA * EffectVarB)
  varcov[3,1] <- HiddenEffectscorrelation * sqrt(EffectVarA * EffectVarB)
  
  varcov[1,2] <- effectSlopecorrelationA * sqrt(EffectVarA * slopeVarA)
  varcov[2,1] <- effectSlopecorrelationA * sqrt(EffectVarA * slopeVarA)
  
  varcov[3,4] <- effectSlopecorrelationB * sqrt(EffectVarB * slopeVarB)
  varcov[4,3] <- effectSlopecorrelationB * sqrt(EffectVarB * slopeVarB)
  
  varcov[2,4] <- slopescorrelation * sqrt(slopeVarA * slopeVarB)
  varcov[4,2] <- slopescorrelation * sqrt(slopeVarA * slopeVarB)
  
  data <- Datagen(n_samples, sample_size, means=c(HiddenEffectA, slopeA, HiddenEffectB, slopeB), varcov = varcov, seed)
  
  saveRDS(data, file = paste(dirpath,"/",HiddenEffectA,"_",EffectVarA, "_", slopeA, "_", slopeVarA,"_", HiddenEffectB,"_",EffectVarB, "_", slopeB, "_", slopeVarB, "_", HiddenEffectscorrelation, "_", slopescorrelation, "_", effectSlopecorrelationA, "_", effectSlopecorrelationB, ".rds", sep = ""))
}

# declare the ranges of data.
HiddenEffect <- c(1,2,3,4)
EffectVar <- c(1,3,5)
slope <-seq(0.8,1.1,1.7)
slopeVar <- c(0.2^2,0.7^2)
HiddentEffectscorrelation <- c(0.2, 0.35)
slopescorrelation <- c(0, 0.1)
effectSlopecorrelation <- c(-0.35, 0, 0.35)

# Create all possible combinations of variables
varCombinations <- expand.grid(HiddenEffect, HiddenEffect, EffectVar, EffectVar, slope, slope, slopeVar, slopeVar, HiddentEffectscorrelation,slopescorrelation,effectSlopecorrelation,effectSlopecorrelation)

# If we're running with MPI, apply in parallel, otherwise apply in serial
if(run_mpi){
  success <- mpi.parApply(varCombinations, 1, gen_data_file, dirpath=path, n_samples=n_samples)
  mpi.close.Rslaves() 
  mpi.quit()
} else {
  success <- apply(varCombinations, 1, gen_data_file, dirpath=path, n_samples=n_samples)
}

