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
    mpi.setup.rngstream(12345)
    mpi.bcast.Robj2slave(all=TRUE)
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
Datagen <- function(n_samples, sample_size, groups,mean,vars, seed){
  
  #create a sample_size x n_samples array x groups
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  sample <- array(rmnorm(n = n_samples*sample_size,mean = mean, varcov = diag(vars)), c(sample_size, n_samples, groups))
  
  #now change it to a more intelligible sample_size x groups x n_samples array
  sample <- aperm(sample, c(1,3,2))
  
  sample[,1,] <- sample[,1,]*sample[,2,]
  
  return(sample)
  
}

# create data and save to file function
gen_data_file <- function(varVector, dirpath, n_samples=1000, sample_size=30, seed = 1234){
  
  HiddenEffect <- varVector[1]
  EffectVar <- varVector[2]
  slope <- varVector[3]
  slopeVar <- varVector[4]
  
  data <- Datagen(n_samples, sample_size, groups = 2, mean=c(HiddenEffect, slope), vars = c(EffectVar, slopeVar), seed)
  
  saveRDS(data, file = paste(dirpath,"/",HiddenEffect,"_",EffectVar, "_", slope, "_", slopeVar, ".rds", sep = ""))
}

# declare the ranges of data.
HiddenEffect <- c(1,2,3)
EffectVar <- c(0.1,0.5)
slope <-seq(0.6,1.4)
slopeVar <- c(0.01,0.05)

# Create all possible combinations of variables
varCombinations <- expand.grid(HiddenEffect, EffectVar, slope, slopeVar)

# If we're running with MPI, apply in parallel, otherwise apply in serial
if(run_mpi){
  success <- mpi.parApply(varCombinations, 1, gen_data_file, dirpath=path)
} else {
  success <- apply(varCombinations, 1, gen_data_file, dirpath=path)
}

