## ---------------------------------- 
## wisdm - parameter tuning functions  
## ApexRMS, May 2025       
## ---------------------------------- 

library(tidyverse)
library(png)
if(!"magick" %in% rownames(installed.packages())){install.packages("magick")}
library(magick)

# Functions --------------------------------------------------------------------

## Build Placeholder Image Function ----

# example inputs
# outputPath <- file.path(out$tempDir, "UnableToFitModel.png")
  
buildFillerImage <- function(outputPath){
 # build replacement png for model output files
  png(outputPath, width = 1000, height = 1000, bg = "white")  # or bg = "transparent"
  par(mar = c(0, 0, 0, 0))  # Remove margins
  plot.new()                # Start new plot
  plot.window(xlim = c(0, 1), ylim = c(0, 1))  # Set coordinate system
  text(0.5, 0.5, "Unable to fit model with defined parameters.", cex = 2, col = "black", font = 2)
  dev.off() 

  } # end funciton


## Build Tuning Matrix Function ----

# This script generates matrices of the wisdm image outputs, where each dimension
# displays a parameter that varied across scenarios.

# example inputs 
# modType <- modType
# modelOutputsTable <- modelOutputsTable
# parameters <- parameterNames
# outputPath <- ssimTempDir  

buildTuningMatrices <- function(modType, 
                                modelOutputsTable, 
                                parameters, 
                                outputPath){
  
  # list of model outputs that matrices can be built for
  modelOutputs <- c("ResponseCurves",
                    "ResidualsPlot",
                    "ResidualsSmoothPlot",
                    "CalibrationPlot",
                    "ROCAUCPlot",
                    "AUCPRPlot",
                    "ConfusionMatrix",
                    "VariableImportancePlot")
  
  
  # Crosswalk of parameter names to clean names
  if(modType == "brt"){
    parameterCW <- tibble(
      Name = c("LearningRate", "NumberOfTrees", "BagFraction", "MaximumTrees"),
      CleanName = c("Learning Rate", "Number of Trees Added Per Stage", "Bag Fraction", "Maximum Number of Trees"))
  } else if(modType == "rf"){
    parameterCW <- tibble(
      Name = c("NumberOfVariablesSampled", "MaximumNodes", "NumberOfTrees", "NodeSize"),
      CleanName = c("Number of Variables Sampled at Split", "Maximum Number of Nodes", "Number of Trees", "Node Size"))
  }
  
  # Get model outputs
  modelOutputsTable <- modelOutputsTable %>%
    dplyr::select(-"ModelsID", -"ModelRDS", -"ResidualSmoothRDS", -"TextOutput", -"VariableImportanceData", -"MaxentFiles") %>%
    pivot_longer(cols = seq(length(parameters)+1,ncol(.)),
                 names_to = "ModelOutput",
                 values_to = "ImageFile") %>%
    right_join(data.frame(ModelOutput = modelOutputs))
  
  # Create output matrices
  for(modelOutput in modelOutputs){
    
    # Get data for model output 
    modelOutputData <- modelOutputsTable %>% 
      filter(ModelOutput == modelOutput)
    
    # Check if data is available for the modelOutput type
    if(modelOutputData$ImageFile %>% unique() %>% length == 1){
      if(is.na(modelOutputData$ImageFile %>% unique())){
        print(str_c(modelOutput, " was not a wisdm output."))
      }
    } else {
      
      # Get # of columns and rows for matrix
      # - NB: if only one parameter varied across scenarios, outputs are arranged in columns
      if(length(parameters) > 1){
        
        numberOfColumns <- modelOutputData %>% 
          dplyr::select(parameters[2] %>% all_of()) %>% 
          n_distinct()
        
        columnNames <- modelOutputData %>% 
          pull(parameters[2] %>% all_of()) %>% 
          unique() %>% 
          as.character()
        
        numberOfRows <- modelOutputData %>% 
          dplyr::select(parameters[1] %>% all_of()) %>% 
          n_distinct()
        
        rowNames <- modelOutputData %>% 
          pull(parameters[1] %>% all_of()) %>% 
          unique() %>% 
          as.character()
        
      } else {
        numberOfColumns <- modelOutputData %>% 
          dplyr::select(parameters[1] %>% all_of()) %>% 
          n_distinct() 
        
        columnNames <- modelOutputData %>% 
          pull(parameters[1] %>% all_of()) %>% 
          unique() %>% 
          as.character()
      }
      
      # Read wisdm png files
      imageList <- image_read(modelOutputData$ImageFile)
      
      # Get width and height of an image
      imageWidth <- imageList[1] %>% 
        image_info() %>% 
        pull(width)
      
      imageHeight <- imageList[1] %>% 
        image_info() %>% 
        pull(height)
      
      # Define font sizes of annotations
      individualPanelFontSize <- imageHeight * 0.04
      matrixFontSize <- imageHeight * 0.02
      
      # Annotate images with column and row names
      # - NB: if only one parameter varied across scenarios, output are arrange in columns
      if(length(parameters) > 1){
        
        # Create a white column to add space for row annotations
        whiteColumn <- image_blank(width = 200,
                                   height = imageHeight,
                                   color = "white")
        
        # Add white space to the left edge of images in first column
        for(i in seq(1:(numberOfRows - 1))){
          
          columnIndex <- (i * numberOfColumns) + 1
          
          if(i == 1){
            # Append and annotate the first image
            imageList[i] <- image_append(c(whiteColumn, imageList[i])) %>% 
              image_annotate(text = rowNames[i],
                             gravity = "West",
                             location = "+100+0",
                             degrees = 270,
                             size = individualPanelFontSize)
            
            # Append and annotate the first image in the second row
            imageList[columnIndex] <- image_append(c(whiteColumn, imageList[columnIndex])) %>% 
              image_annotate(text = rowNames[i + 1],
                             gravity = "West",
                             location = "+100+0",
                             degrees = 270,
                             size = individualPanelFontSize)
          } else
            # Append the first image of remaining rows
            imageList[columnIndex] <- image_append(c(whiteColumn, imageList[columnIndex])) %>% 
              image_annotate(text = rowNames[i + 1],
                             gravity = "West",
                             location = "+100+0",
                             degrees = 270,
                             size = individualPanelFontSize)
        }
      }
      
      # Create a white row to add space for column annotations
      whiteRow <- image_blank(width = imageWidth,
                              height = 200,
                              color = "white")
      
      # Add white space to the top edge of all images to ensure consistent image sizes
      for(i in seq(1:length(imageList))){
        imageList[i] <- image_append(c(whiteRow, imageList[i]), stack = TRUE)
      }
      
      
      # Annotate images in first row with column names
      imageList[1:numberOfColumns] <- imageList[1:numberOfColumns] %>% 
        image_annotate(text = columnNames,
                       gravity = "North",
                       location = "+0+100",
                       size = individualPanelFontSize)
      
      # Define matrix tiles
      if(length(parameters) > 1){
        maxtrixTiles <- numberOfColumns %>% as.character() %>% 
          str_c("x", numberOfRows %>% as.character())
      } else {
        maxtrixTiles <- numberOfColumns %>% as.character() %>% 
          str_c("x1")
      }
      
      # Define matrix panel sizes and spacing
      matrixGeometry <- "x" %>% 
        str_c((imageWidth/numberOfColumns) %>% as.character(),
              "+0+0")
      
      # Create the montage
      montage_image <- imageList %>%
        image_montage(geometry = matrixGeometry, tile = maxtrixTiles)
      montage_image_info <- image_info(montage_image)
      
      if(length(parameters) > 1){ 
        
        # Create a blank image with extra space 
        blank_image <- image_blank(width = montage_image_info$width+50, height = montage_image_info$height+50, color = "white")
        
        # Create matrix
        outputMatrix <- blank_image %>%
          image_composite(montage_image, offset = "+50+50") %>% 
        # outputMatrix <- imageList %>% 
        #   image_montage(geometry = matrixGeometry, tile = maxtrixTiles) %>% 
          image_annotate(
            text = parameterCW %>% filter(Name == parameters[2]) %>% pull(CleanName),
            gravity = "North",
            location = "+0+25",
            size = matrixFontSize)%>% 
          image_annotate(
            text = parameterCW %>% filter(Name == parameters[1]) %>% pull(CleanName),
            gravity = "West",
            location = "+25+200",
            degrees = 270,
            size = matrixFontSize)
      } else {
        # Create matrix
        outputMatrix <- imageList %>% 
          image_montage(geometry = matrixGeometry, tile = maxtrixTiles) %>% 
          image_annotate(
            text = parameterCW %>% filter(Name == parameters[1]) %>% pull(CleanName),
            gravity = "North",
            location = "+0+0",
            size = matrixFontSize)
      }
      
      outputName <- # modType %>% str_to_lower() %>% 
        str_c(modelOutput, "Matrix.png")
      
      if(!file.exists(file.path(outputPath, outputName))){
        image_write(outputMatrix,
                    path = file.path(outputPath, outputName),
                    format = "png",
                    quality = 100)
      } else stop("The output matrix already exists. Delete the file before rewriting to disk.")
    }
  }
} # end function

## Combine text files function ----

  combineTxtFiles <- function(filePaths,  # list of file paths for .txt files to be combined  
                              outputPath  # file path for .txt output
  ){
    # Read and combine all files
    combined_text <- sapply(filePaths, function(x) paste(readLines(x), collapse = "\n"))
    
    # Combine all into one big string
    final_text <- paste(combined_text, collapse = "\n\n============================================================\n\n")
    
    # Write to a new file
    writeLines(final_text, outputPath)
    
    } # end function 