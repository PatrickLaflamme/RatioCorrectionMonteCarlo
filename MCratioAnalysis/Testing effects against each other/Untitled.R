#library(Rmpi)
library(mnormt)

Datagen <- function(n_samples, sample_size, groups,mean,vars, seed){
  
  #create a sample_size x n_samples array x groups
  
  if(!is.null(seed)){
    set.seed(seed)
  }
  sample <- array(rmnorm(n = n_samples*sample_size,mean = mean, varcov = diag(vars)), c(sample_size, n_samples, groups))
  
  #now change it to a more intelligible sample_size x groups x n_samples array
  sample <- aperm(sample, c(1,3,2))
  
  return(sample)
  
}

gen_data_file <- function(HiddenEffect, EffectVar, slope, slopeVar, dirpath, n_samples=100000, sample_size=30, seed = 1234){
  
  data <- Datagen(n_samples, sample_size, groups = 2, mean=c(HiddenEffect, slope), vars = c(EffectVar, slopeVar), seed)
  
  out <- system.time(save(file = paste(dirpath,"/",HiddenEffect,"_",EffectVar, "_", slope, "_", slopeVar, ".Rdata", sep = ""), "data"))
  
  print(out)
}

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

HiddenEffect <- seq(0,10,by=1)

EffectVar <- c(0.1,0.5,1,1.5,2)

slope <-seq(0.6,1.4, by=0.1)

slopeVar <- c(0.01,0.05,0.1,0.2,0.3)

varCombinations <- expand.grid(HiddenEffect, EffectVar, slope, slopeVar)

colnames(varCombinations) <- c("Hidden Effect", "Effect variance", "Slope", "slope Variance")

