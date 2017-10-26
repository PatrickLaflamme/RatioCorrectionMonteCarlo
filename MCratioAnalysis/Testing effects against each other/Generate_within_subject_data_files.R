#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#check what the commandArgs are requesting
if(length(args) > 0){
  
  # Should we run with MPI?
  run_mpi <- length(grep("mpi", args))
  
  # Where to save the data?
  path <- grep("--path", args) + 1
  
  # Where to save the data?
  n_samples <- grep("--n_samples", args) + 1
  
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
  if( length(n_samples) ) {
    n_samples <- args[n_samples]
  } else {
    stop("n_samples must be provided using --n_samples [integer]")
  }
} else {
  
  run_mpi <- 0
  
  path <- '~/Desktop/Test/'
}

library(mnormt)

# data generation function
Datagen <- function(n_samples, sample_size, means,varcov, seed){
  
  #create a sample_size x n_samples array x groups
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  sample <- array(rmnorm(n = n_samples*sample_size,mean = means, varcov = varcov), c(sample_size, n_samples, length(means)))
  
  #now change it to a more intelligible sample_size x groups x n_samples array
  sample <- aperm(sample, c(1,3,2))
  
  return(sample)
  
}

# create data and save to file function
gen_data_file <- function(varVector, dirpath, n_samples=n_samples, sample_size=30, seed = 1234){
  
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
  
  varcov <- diag(c(EffectVarA, EffectVarB,slopeVarA,slopeVarB))
  
  varcov[1,2] <- HiddenEffectscorrelation * sqrt(EffectVarA * EffectVarB)
  varcov[2,1] <- HiddenEffectscorrelation * sqrt(EffectVarA * EffectVarB)
  
  varcov[1,3] <- effectSlopecorrelationA * sqrt(EffectVarA * slopeVarA)
  varcov[3,1] <- effectSlopecorrelationA * sqrt(EffectVarA * slopeVarA)
  
  varcov[2,4] <- effectSlopecorrelationB * sqrt(EffectVarB * slopeVarB)
  varcov[4,2] <- effectSlopecorrelationB * sqrt(EffectVarB * slopeVarB)
  
  varcov[3,4] <- slopescorrelation * sqrt(slopeVarA * slopeVarB)
  varcov[4,3] <- slopescorrelation * sqrt(slopeVarA * slopeVarB)
  
  data <- Datagen(n_samples, sample_size, means=c(HiddenEffectA, slopeA, HiddenEffectB, slopeB), varcov = varcov, seed)
  
  saveRDS(data, file = paste(dirpath,"/",HiddenEffectA,"_",EffectVarA, "_", slopeA, "_", slopeVarA,"_", HiddenEffectB,"_",EffectVarB, "_", slopeB, "_", slopeVarB, "_", HiddenEffectscorrelation, "_", slopescorrelation, "_", effectSlopecorrelationA, "_", effectSlopecorrelationB, ".rds", sep = ""))
}

# declare the ranges of data.
HiddenEffect <- c(1,2,3,4)
EffectVar <- c(1,3,5)
slope <-seq(0.8,1.1,1.7)
slopeVar <- c(0.1^2,0.7^2)
HiddentEffectscorrelation <- c(0.2, 0.35)
slopescorrelation <- c(0, 0.1)
effectSlopecorrelation <- c(-0.35, 0, 0.35)

# Create all possible combinations of variables
varCombinations <- expand.grid(HiddenEffect, HiddenEffect, EffectVar, EffectVar, slope, slope, slopeVar, slopeVar, HiddentEffectscorrelation,slopescorrelation,effectSlopecorrelation,effectSlopecorrelation)

# If we're running with MPI, apply in parallel, otherwise apply in serial
if(run_mpi){
  success <- mpi.parApply(varCombinations, 1, gen_data_file, dirpath=path)
  mpi.close.Rslaves() 
  mpi.quit()
} else {
  success <- apply(varCombinations, 1, gen_data_file, dirpath=path)
}

