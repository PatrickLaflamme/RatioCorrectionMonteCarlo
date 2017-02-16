skewness <- function (x, na.rm = FALSE, type = 3) 
{
    if (any(ina <- is.na(x))) {
        if (na.rm) 
            x <- x[!ina]
        else return(NA)
    }
    if (!(type %in% (1:3))) 
        stop("Invalid 'type' argument.")
    n <- length(x)
    x <- x - mean(x)
    y <- sqrt(n) * sum(x^3)/(sum(x^2)^(3/2))
    if (type == 2) {
        if (n < 3) 
            stop("Need at least 3 complete observations.")
        y <- y * sqrt(n * (n - 1))/(n - 2)
    }
    else if (type == 3) 
        y <- y * ((1 - 1/n))^(3/2)
    y
}

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
    contour(x=x, y=y, z=contourz, add=T, levels=c(0.045,0.055), col = 'white', lwd=1.5)
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

Datagen <- function(n_samples, sample_size, groups,mean,varcov, seed){
  
  library(mnormt)
  #create a sample_size x n_samples array x groups
  
  if(!is.null(seed)){
  set.seed(seed)
  }
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
  return(indexData)
}

runTest <- function(slopemean, slopevar, effectvar, r, seed = NULL,  sample_size=20, alpha=0.05, n_samples=100000, plot=F, test=T){
  ## n_samples    :: Number of times to simulate the experiment
  ##
  ## sample_size  :: Number of participants in each condition of the experiment
  ##
  ## varcov       :: Variance-covariance matrix of the bivariate normal data source
  ##
  ## mean         :: Population means of the two variables
  ##
  ## condition    :: list of conditions to generate
  ##
  ##                 type 1: participant-wise difference/slope
  ##                 type 2: participant-wise difference/overall slope
  ##                 type 3: taylor approximation of the fieler correction
  ##
  ## alpha        :: The nominal type 1 error rate.
  ##
  
  covar <- r*sqrt(slopevar)*sqrt(effectvar)
  
  varcov <- array(c(effectvar, covar, covar, slopevar), dim=c(2,2))
  
  data <- Datagen(n_samples, sample_size, groups = 2, c(0, slopemean), varcov, seed)
  
  popIndex <- 0
  
  indexed <- list()
  
  for(i in 1:3){
    indexed <- c(indexed, list(genIndex(data, i) - popIndex))
  }
  
  if(test){
    results <- c()
    for(i in 1:3){
    results <- cbind(results, apply(indexed[[i]], 2, mean))
    }
    
    return(results)
  }
  
  return(indexed)
  
}
