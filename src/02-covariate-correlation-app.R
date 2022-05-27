ui <- fluidPage(
  
  # App title ----
  titlePanel("Covariate Correlation Viewer"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(width = 3, 
                 
                 # Input: select which variables to consider 
                 checkboxGroupInput(inputId = "show_vars", 
                                    label = "Covariates to include:",
                                    choices = names(covData), 
                                    selected = names(covData)),
                 # Input: Numeric entry for minimum correlation ----
                 numericInput(inputId = "minCor",
                              label = "Correlation Threshold:",
                              value = options$CorrelationThreshold,
                              min = 0, max = 1, step = 0.05),
                 
                 # Input: Numeric entry for number of plots to view ----
                 numericInput(inputId = "numPlots",
                              label = "Number of Plots:",
                              value = options$NumberOfPlots, 
                              min = 1, max = ncol(siteData)-1, step = 1),
                 
                 # action button
                 actionButton(inputId = "update", label = "Update")
                 
    ),
    
    # Main panel for displaying outputs ----
    mainPanel( width = 9,
               # Output: Covariate correlation image
               imageOutput(outputId = "image", height = "auto", width = "auto")
    )
  )
)


server <- function(input, output, session) { # 
  
  # create a temp file to save the output
  outfile <- file.path(ssimTempDir, "CovariateCorrelationMatrix.png")
  
  v <- reactiveValues(options = options, 
                      selectedCovs = names(covData))
  
  # image creates a new PNG file each time inputs change
  observeEvent(input$update, {
    v$options$CorrelationThreshold <- input$minCor
    v$options$NumberOfPlots <- input$numPlots
    SelectedCovariates <<- v$selectedCovs <- input$show_vars
  })
  
  output$image <- renderImage({
    
    # options$CorrelationThreshold  <-  
    # options$NumberOfPlots <- input$numPlots
    # selectedCovs <- input$show_vars
    
    # Generate the image and write it to file
    pairsExplore(inputData = siteData,
                 selectedCovs = v$selectedCovs,
                 options = v$options,
                 factorVars = factorVars,
                 family = modelFamily,
                 outputFile = outfile)
    
    # Return a list containing information about the image
    list(src = outfile,
         contentType = "image/png",
         height = "100%",
         width = "100%",
         alt = "This is alternate text")
  }, deleteFile = FALSE)
  
  session$onSessionEnded(function() {
    stopApp()
  })
}

shinyApp(ui, server)