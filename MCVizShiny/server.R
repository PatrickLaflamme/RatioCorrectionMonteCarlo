
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

load('MCresults-Labelled.RData')
source('Monte Carlo Basic Functions.R')


fselect_in <- function(x, ref, d = 10){
  round(x, digits=d) %in% round(ref, digits=d)
}

shinyServer(function(input, output) {
  
  output$distPlot <- renderPlot({
    
    if(input$datatype=='Skewness-Pop'){
      load('MCresults-Skewness.RData')
      dimnames(output) <- dimnames(newOutput)
      newOutput <- output
      ERtick <- pretty(newOutput, n=5)
      cont <- NA
    }
    if(input$datatype=='Skewness-SampDist'){
      load('MCresults-Skewness-SampDist.RData')
      dimnames(output) <- dimnames(newOutput)
      newOutput <- output
      ERtick <- pretty(newOutput, n=5)
      cont <- NA
    }
    else if(input$datatype=='Actual Alpha Values'|input$datatype=='Alphas and SampDist Skew'){
      load('MCresults-Labelled.RData')
      ERtick <- c(0.005, 0.01,0.025,0.05,0.1,0.2,0.4,0.8)
    }
    
    requestedcor <- fselect_in(as.numeric(dimnames(newOutput)[[2]]),input$correlation, d=2)
    requestedmean <- fselect_in(as.numeric(dimnames(newOutput)[[3]]),input$slopemean, d=2)
    
    if(all(ERtick == c(0.005, 0.01,0.025,0.05,0.1,0.2,0.4,0.8))){
      cont <- cont <- unname(newOutput[input$method,dimnames(newOutput)[[2]][requestedcor],dimnames(newOutput)[[3]][requestedmean],,])
    }
    
    MCfilled.contour(x=seq(from=0,to=10,length.out=20),
                     y=seq(from=0,to=1,length.out=20),
                     z=unname(newOutput[input$method,dimnames(newOutput)[[2]][requestedcor],dimnames(newOutput)[[3]][requestedmean],,]),
                     zlim = c(min(newOutput),max(newOutput)),
                     contourz = cont,
                     color.palette = MCrainbow,
                     nlevels = 100,
                     xlab = 'Effect Variance',
                     ylab = 'Slope Variance',
                     main = paste(input$method,input$correlation,input$slopemean,'alpha = 0.05', sep=' '),
                     key.axes = axis(4, at = c(ERtick, round(max(newOutput), 3), label = c(ERtick, round(max(newOutput), 3)))),
                     key.title = expression(paste('Actual ', alpha, sep = '')))

  })
  
  output$secPlot <- renderPlot({
    
    if(input$datatype=='Alphas and SampDist Skew'){
      load('MCresults-Skewness-SampDist.RData')
      dimnames(output) <- dimnames(newOutput)
      skewtick <- pretty(output, n=5)
    }
    else{return(NULL)}
    
    requestedcor <- fselect_in(as.numeric(dimnames(output)[[2]]),input$correlation, d=2)
    requestedmean <- fselect_in(as.numeric(dimnames(output)[[3]]),input$slopemean, d=2)
    
    MCfilled.contour(x=seq(from=0,to=10,length.out=20),
                     y=seq(from=0,to=1,length.out=20),
                     z=unname(output[input$method,dimnames(output)[[2]][requestedcor],dimnames(output)[[3]][requestedmean],,]),
                     zlim = c(min(output),max(output)),
                     color.palette = MCrainbow,
                     nlevels = 100,
                     xlab = 'Effect Variance',
                     ylab = 'Slope Variance',
                     main = paste(input$method,input$correlation,input$slopemean,'alpha = 0.05', sep=' '),
                     key.axes = axis(4, at = c(skewtick, round(max(output), 3), label = c(skewtick, round(max(output), 3)))),
                     key.title = expression(paste('Actual ', alpha, sep = '')))
  })

})
