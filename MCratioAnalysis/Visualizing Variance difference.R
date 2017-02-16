

val0 <- runSkew(5,0.9, 0.2, 1, 0, test=F)
#test0 <- runTest(5,0.9, 0.2, 1, 0)

par(mfrow=c(2,1), oma=c(0, 0, 2, 0))
hist(apply(val0[[2]],2,var),breaks = seq(0,1000,by = 0.1), xlim = c(0,50),main = "ZeroVar", xlab = 'Variance')
hist(apply(val0[[3]],2,var),breaks = seq(0,1000,by = 0.1), xlim = c(0,50),main = "Taylor", xlab = 'Variance')
title(main = "Histograms: E=5, S=0.9, varE=1, varS=0.2, r=0",outer=T)

val0.4 <- runSkew(5,0.9, 0.2, 1, 0.4, test=F)
#test0.4 <- runTest(5,0.9, 0.2, 1, 0.4)

par(mfrow=c(2,1), oma=c(0, 0, 2, 0))
hist(apply(val0.4[[2]],2,var),breaks = seq(0,1000,by = 0.1), xlim = c(0,50),main = "ZeroVar", xlab = 'Variance')
hist(apply(val0.4[[3]],2,var),breaks = seq(0,1000,by = 0.1), xlim = c(0,50),main = "Taylor", xlab = 'Variance')
title(main = "Histograms: E=5, S=0.9, varE=1, varS=0.2, r=0.4",outer=T)

