#library(R.matlab)

#writeMat('MCresults.mat', MCresults = output)

system('mkdir MCresultPlots')


effectvarRange <- seq(from=0, to=10, length.out = 20)
slopevarRange <- seq(from=0, to=1, length.out = 20)
ERtick <- c(0.005, 0.01,0.025,0.05,0.1,0.2,0.4,0.8)



for(i in dimnames(output)[[1]]){
  
  system(paste('mkdir MCresultPlots/', i, sep=''))
  
  for(j in dimnames(output)[[2]]){
    
    system(paste('mkdir MCresultPlots/', i, '/', j, sep=''))
    
    for(k in dimnames(output)[[3]]){
      
      folder <- paste('MCresultPlots/', i, '/', j, '/',  sep='')
      
      pdf(file=paste(folder, k,'.pdf',sep=''), width = 10, height=10)
      MCfilled.contour(x=effectvarRange,
                       y=slopevarRange,
                       z=unname(output[i,j,k,,]),
                       zlim = c(0,0.1),
                       contourz = unname(output[i,j,k,,]),
                       color.palette = MCrainbow,
                       nlevels = 100,
                       xlab = 'Effect Variance',
                       ylab = 'Slope Variance',
                       main = paste(i,j,k,expression(alpha),'= 0.05', sep=' '),
                       key.axes = axis(4, at = ERtick, label = ERtick),
                       key.title = expression(paste('Actual ', alpha, sep = ''))
                      )
      dev.off()
    }
    
  }
  
}
