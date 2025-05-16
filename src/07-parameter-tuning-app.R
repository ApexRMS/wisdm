## ----------------------------------
## wisdm - parameter tuning shiny app
## ApexRMS, May 2025
## ----------------------------------


ui <- fluidPage(
  
  # App title ----
  titlePanel("Parameter Tuning Viewer"),
  
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    
    # Sidebar panel for inputs ----
    sidebarPanel(width = 3, 
                 
                 # Input: select which parameters to consider 
                 radioButtons(inputId = "selectedCombo", 
                              label = "Parameter combinations:",
                              choices = displayNames, 
                              selected = NULL), 
                 # # Input: Numeric entry for number of plots to view ----
                 # numericInput(inputId = "numPlots",
                 #              label = "Number of Plots:",
                 #              value = nrow(comboImgs), # options$NumberOfPlots, 
                 #              min = 1, max = nrow(comboImgs), step = 1),
                 # 
                 # # action button
                 # # actionButton(inputId = "update", label = "Update"),
                 
                 # action button
                 actionButton(inputId = "close", label = "Save & Close"),
                 tags$script( "Shiny.addCustomMessageHandler('closeWindow', 
                              function(data) {eval(data.message)});" )
    ),
    
    # Main panel for displaying outputs ----
    mainPanel(width = 9,
              tabsetPanel(
                tabPanel("Response Curves", imageOutput(outputId = "image1", height = "auto", width = "auto")),
                tabPanel("Calibration", imageOutput(outputId = "image2", height = "auto", width = "auto")),
                tabPanel("ROCAUC", imageOutput(outputId = "image3", height = "auto", width = "auto")),
                tabPanel("AUCPR", imageOutput(outputId = "image4", height = "auto", width = "auto")),
                tabPanel("Confusion Matrix", imageOutput(outputId = "image5", height = "auto", width = "auto")),
                tabPanel("Variable Importance", imageOutput(outputId = "image6", height = "auto", width = "auto")),
                ) 
              )
    )
)


server <- function(input, output, session) { # 
  
  v <- reactiveValues(selectedCombo = displayNames) 
  
  # # image creates a new PNG file each time inputs change
  # observeEvent(input$update, {
  #   v$options$NumberOfPlots <- input$numPlots
  #   v$selectedCombo <- input$selectedCombo 
  # })
  
  output$image1 <- renderImage({
    list(src = file.path(ssimTempDir, "ResponseCurvesMatrix.png"), width = "100%")
    }, deleteFile = FALSE)
  
  output$image2 <- renderImage({
    list(src = file.path(ssimTempDir, "CalibrationPlotMatrix.png"), width = "100%")
  }, deleteFile = FALSE)
  
  output$image3 <- renderImage({
    list(src = file.path(ssimTempDir, "ROCAUCPlotMatrix.png"), width = "100%")
  }, deleteFile = FALSE)
  
  output$image4 <- renderImage({
    list(src = file.path(ssimTempDir, "AUCPRPlotMatrix.png"), width = "100%")
  }, deleteFile = FALSE)
  
  output$image5 <- renderImage({
    list(src = file.path(ssimTempDir, "ConfusionMatrixMatrix.png"), width = "100%")
  }, deleteFile = FALSE)
  
  output$image6 <- renderImage({
    list(src = file.path(ssimTempDir, "VariableImportancePlotMatrix.png"), width = "100%")
  }, deleteFile = FALSE)
  
  observeEvent(input$close, {
    comboOut <<- input$selectedCombo 
    session$sendCustomMessage(type = "closeWindow", list(message ="window.close()"))
    stopApp()
  })
  
  session$onSessionEnded(function() {
    stopApp()
  })
}

shinyApp(ui, server)