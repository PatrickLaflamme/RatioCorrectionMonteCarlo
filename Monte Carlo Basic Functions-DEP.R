MCfilled.contour <- function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                                           length.out = ncol(z)), z, contourz=NULL, xlim = range(x, finite = TRUE), 
                              ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
                              levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
                              col = color.palette(length(levels) - 1), plot.title, plot.axes, 
                              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
                              axes = TRUE, frame.plot = axes, border=NA, ...) 
{
  if (missing(z)) {
    if (!missing(x)) {
      if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
      }
      else {
        z <- x
        x <- seq.int(0, 1, length.out = nrow(z))
      }
    }
    else stop("no 'z' matrix specified")
  }
  else if (is.list(x)) {
    y <- x$y
    x <- x$x
  }
  if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
    stop("increasing 'x' and 'y' values expected")
  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
  on.exit(par(par.orig))
  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
  layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
  par(las = las)
  mar <- mar.orig
  mar[4L] <- mar[2L]
  mar[2L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
              yaxs = "i")
  rect(0, levels[-length(levels)], 1, levels[-1L], col = col, border = border)
  if (missing(key.axes)) {
    if (axes) 
      axis(4)
  }
  else key.axes
  box()
  if (!missing(key.title)) 
    key.title
  mar <- mar.orig
  mar[4L] <- 1
  par(mar = mar)
  plot.new()
  plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
  .filled.contour(x, y, z, levels, col)
  if(!missing(contourz))
    contour(x=x, y=y, z=contourz, add=T, levels=c(0.04,0.06), col = 'white', lwd=1.5)
  if (missing(plot.axes)) {
    if (axes) {
      title(main = "", xlab = "", ylab = "")
      Axis(x, side = 1)
      Axis(y, side = 2)
    }
  }
  else plot.axes
  if (frame.plot) 
    box()
  if (missing(plot.title)) 
    title(...)
  else plot.title
  invisible()
}

MCrainbow <-function (n, s = 1, v = 1, start = 0.36, end = 0.15, 
                      alpha = 1) 
{
  if ((n <- as.integer(n[1L])) > 0) {
    if (start == end || any(c(start, end) < 0) || any(c(start, 
                                                        end) > 1)) 
      stop("'start' and 'end' must be distinct and in [0, 1].")
    hsv(h = seq.int(start, ifelse(start > end, 1, 0) + end, 
                    length.out = n)%%1, s, v, alpha)
  }
  else character()
}

Datagen <- function(n_samples, sample_size, groups,mean,varcov){
  
  library(mnormt)
  #create a sample_size x n_samples array x groups
  sample <- array(rmnorm(n = n_samples*sample_size,mean = mean, varcov = varcov), c(sample_size, n_samples, groups))
  
  #now change it to a more intelligible sample_size x groups x n_samples array
  sample <- aperm(sample, c(1,3,2))
  
  return(sample)
  
}

genIndex <- function(data, condition){
  ## Generates an index score for each data point. Generates it in one of 3 ways, as listed in the 'type' of condition.
  ##
  ## data      :: always has shape: sampleSize x 2 x n_samples
  ##                                The 2 columns are: random sample dY and slope
  ## condition :: type 1: participant-wise difference/slope
  ##              type 2: participant-wise difference/overall slope
  ##              type 3: taylor-corrected effect scores based on Feiller's theorem
  
  assertthat::are_equal(dim(data)[2], 2)
  
  if(condition == 1){
  indexData <- data[,1,]/data[,2,] #divide each individual effect by each individual slope
  }
  
  else if(condition ==2){
    indexData <- t(t(data[,1,])/apply(data[,2,], 2, mean)) #divide each individual effect by the sample mean slope
  }
  
  else if (condition==3){
    meanSlope <- apply(data[,2,], 2, mean) #calculate the sample mean slope
    meanEffect <- apply(data[,1,], 2, mean) #calculate the sample mean effect
    
    indexData <- meanEffect/meanSlope * (1 + t(t(data[,1,])/meanEffect) - t(t(data[,2,])/meanSlope)) #perform the correction
  }
  return(indexData)
}

runTest <- function(n_samples, sample_size, varcov, mean, condition, alpha, plot=F){
  ## n_samples    :: Number of times to simulate the experiment
  ##
  ## sample_size  :: Number of participants in each condition of the experiment
  ##
  ## varcov       :: Variance-covariance matrix of the bivariate normal data source
  ##
  ## mean         :: Population means of the two variables
  ##
  ## condition    :: type 1: participant-wise difference/slope
  ##                 type 2: participant-wise difference/overall slope
  ##
  ## alpha        :: The nominal type 1 error rate.
  ##
  
  data <- Datagen(n_samples, sample_size, groups = 2, mean, varcov)
  
  popIndex <- mean[1]/mean[2]
  
  indexed <- genIndex(data, condition) - popIndex
  
  if(plot){hist(indexed, density = T)}
  
  results <- apply(indexed, 2, t.test, conf.level = 1 - alpha)
  
  pvalues <- unlist(sapply(results,'[',3))
  
  type1ER <- sum(pvalues<=alpha)/n_samples
  
  return(type1ER)
  
}



# ER1 <- array(data = NA, dim = c(20,20))
# ER2 <- array(data = NA, dim = c(20,20))
# ER3 <- array(data = NA, dim = c(20,20))
# 
# 
# for(i in seq(0.2,4, by = 0.2)){
# 
#   for(j in seq(0.2,4, by = 0.2)){
# 
#     ER1[i*5,j*5] <- runTest(100000, 20, array(c(1,0,0,1), dim=c(2,2)), mean = c(i,j), condition = 1, 0.05)
#     ER2[i*5,j*5] <- runTest(100000, 20, array(c(1,0,0,1), dim=c(2,2)), mean = c(i,j), condition = 2, 0.05)
#     ER3[i*5,j*5] <- runTest(100000, 20, array(c(1,0,0,1), dim=c(2,2)), mean = c(i,j), condition = 3, 0.05)
#     
#     print(c(i*5,j*5))
#   }
# 
# }
# remove(i,j)
# 
# logER1 <- log(ER1)
# logER2 <- log(ER2)
# logER3 <- log(ER3)
# ERtick <- c(0.005, 0.01,0.025,0.05,0.1,0.2,0.4,0.8)
# 
# 
# MCfilled.contour(x = seq(0,3.8,length.out = 20),
#                  y = seq(0.2,4,length.out = 20),
#                  z = logERvar1,
#                  contourz = ERvar1,
#                  zlim = c(-5,0),
#                  nlevels = 100,
#                  color.palette = MCrainbow,
#                  xlab = 'slope variance',
#                  ylab = 'difference variance',
#                  main = expression(paste("Actual ", alpha," for Index method: nominal ",alpha," = 0.05")),
#                  key.axes = axis(4, at = log(ERtick), label = ERtick))
# dev.copy(pdf, '~/Dropbox/MCratio/Index-var.pdf')
# dev.off()
# 
# 
# MCfilled.contour(x = seq(0,3.8,length.out = 20),
#                  y = seq(0.2,4,length.out = 20),
#                  z = logER1,
#                  contourz = ER1,
#                  zlim = c(-5,0),
#                  nlevels = 100,
#                  color.palette = MCrainbow,
#                  xlab = 'slope mean',
#                  ylab = 'effect mean',
#                  main = expression(paste("Actual ", alpha," for Index method: nominal ",alpha," = 0.05")),
#                  key.axes = axis(4, at = log(ERtick), label = ERtick))
# dev.copy(pdf, '~/Dropbox/MCratio/Index-slope-diff.pdf')
# dev.off()
# 
# MCfilled.contour(x = seq(0,3.8,length.out = 20),
#                  y = seq(0.2,4,length.out = 20),
#                  z = logERvar2,
#                  contourz = ERvar2,
#                  zlim = c(-5,0),
#                  nlevels = 100,
#                  color.palette = MCrainbow,
#                  xlab = 'slop variance',
#                  ylab = 'difference variance',
#                  main = expression(paste("Actual ", alpha," for zero-variance method: nominal ",alpha," = 0.05")),
#                  key.axes = axis(4, at = log(ERtick), label = ERtick))
# 
# dev.copy(pdf, '~/Dropbox/MCratio/zeroVariance-var.pdf')
# dev.off()
# 
MCfilled.contour(x = seq(0,3.8,length.out = 20),
                 y = seq(0.2,4,length.out = 20),
                 z = logER3,
                 contourz = ER3,
                 zlim = c(-5.3,0),
                 nlevels = 100,
                 color.palette = MCrainbow,
                 xlab = 'effect mean',
                 ylab = 'slope mean',
                 main = expression(paste("Actual ", alpha," for taylor method: nominal ",alpha," = 0.05")),
                 key.axes = axis(4, at = log(ERtick), label = ERtick))

dev.copy(pdf, '~/Dropbox/MCratio/Taylor-slope-diff.pdf')
dev.off()


