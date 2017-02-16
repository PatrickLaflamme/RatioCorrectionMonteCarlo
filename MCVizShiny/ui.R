
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(

  # Application title
  titlePanel("MC Results Ratio Analysis"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      selectInput('datatype',
                  "Select which data to visualize",
                  choices = c('Alphas and SampDist Skew', 'Skewness-Pop', 'Skewness-SampDist', 'Actual Alpha Values')),
      sliderInput("slopemean",
                  "Select Slope Mean:",
                  min = 0.55,
                  max = 1.45,
                  value = 1.05,
                  step = 0.100000000000001,
                  round=-2),
      sliderInput("correlation",
                  "Select Correlation Value:",
                  min = 0,
                  max = 0.4,
                  value = 1.05,
                  step = (0.4/7),
                  sep = '',
                  round=-16),
      selectInput('method',
                  'Select Analysis Method',
                  choices = c('method-Index', 'method-ZeroVar', 'method-Taylor'))),

    # Show a plot of the generated distribution
    mainPanel(
      plotOutput("distPlot"),
      
      plotOutput('secPlot'),
      
      plotOutput('thirdPlot')
      
    )
  )
))
