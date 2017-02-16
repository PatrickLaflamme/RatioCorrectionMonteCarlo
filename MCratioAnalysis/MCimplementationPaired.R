library(Rmpi)
library(mnormt)

mpi.setup.rngstream(12345)

source('Monte Carlo Basic Functions.R')

mpi.bcast.Robj2slave(all=TRUE)


start <- proc.time()
print(runTest(1,1,1,0,1,1))
print(proc.time()-start)

spaceDims <- c(3,10,10,3,10,3) #3 methods x 10 slopevariances x 10 effect variances x 3 slope means x 10 effect means x 3 correlations

#slope means 0.55 - 1.45 by 0.1
#effect means 0 - 10 by 1
#slope variance 0 - 1 by 0.05
#effect variance 0 - 10 by 0.5
#correlation 0-0.4 by 0.05

starts <- c(0.05,0.5,0.55,1,0)

ends <- c(1,10,1.45,10,0.4)

spacelist <- array( dim = c(prod(spaceDims),length(spaceDims)))

for(i in 2:length(spaceDims)){
  
  curlist <- seq(from = starts[i-1], to = ends[i-1], length.out = spaceDims[i])
  
  spacelist[,i-1] <- rep(curlist, times=prod(spaceDims[-c(1,i:length(spaceDims))]), each = prod(spaceDims[-(1:i)]))
  
}

spacelist <- split(spacelist, seq(nrow(spacelist)))


outputSpaceDims <- c(3,8,10,20,20)

##    spaceDims ::  dimensions of the parameter space that will be explored (range * resolution for each dim)

output <- array(dim = outputSpaceDims)

j<-1
i<-1
print(array(unlist(mpi.apply(spacelist[1600*(i-1)+(j-1)*80+1:80], function(e){ return(runTest(e[3],e[1],e[2], e[4]))})), dim=c(3,8,10)))
print(spacelist[1600*(i-1)+(j-1)*80+1:80])

for(slopevar in 1:20){
  for(effectvar in 1:20){
    output[,,,effectvar,slopevar] <- array(unlist(mpi.apply(spacelist[1600*(slopevar-1)+(effectvar-1)*80+1:80], function(e){ return(runTest(e[3],e[1],e[2], e[4]))})), dim=c(3,8,10))
    
    print(1600*(slopevar-1)+(effectvar)*80)
    
  }
}

save(output, file='MCresults.RData')

mpi.quit(save='yes')
