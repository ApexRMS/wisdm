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

library(tibble)
library(dplyr)
library(tidyr)
library(png)
library(grid)
library(gridExtra)

# buildTuningMatrices: assemble PNG panels into a labeled matrix
buildTuningMatrices <- function(modType,
                                modelOutputsTable,
                                parameters,
                                outputPath) {
  # possible outputs
  modelOutputs <- c(
    "ResponseCurves","ResidualsPlot","ResidualSmoothPlot",
    "CalibrationPlot","ROCAUCPlot","AUCPRPlot",
    "ConfusionMatrix","VariableImportancePlot"
  )
  
  # map raw param names to clean labels
  parameterCW <- if (modType == "brt") {
    tibble(
      Name = c("LearningRate","NumberOfTrees","BagFraction","MaximumTrees"),
      CleanName = c("Learning Rate",
                    "Number of Trees Added Per Stage",
                    "Bag Fraction",
                    "Maximum Number of Trees")
    )
  } else {
    tibble(
      Name = c("NumberOfVariablesSampled","MaximumNodes","NumberOfTrees","NodeSize"),
      CleanName = c("Number of Variables Sampled at Split",
                    "Maximum Number of Nodes",
                    "Number of Trees",
                    "Node Size")
    )
  }
  
  # tidy table: pivot to long, keep only relevant columns
  modelOutputsTable <- modelOutputsTable %>%
    select(-ModelsID, -ModelRDS, -ResidualSmoothRDS,
           -TextOutput, -VariableImportanceData, -MaxentFiles) %>%
    pivot_longer(
      cols = seq(length(parameters)+1, ncol(.)),
      names_to = "ModelOutput", values_to = "ImageFile"
    ) %>%
    right_join(tibble(ModelOutput = modelOutputs), by = "ModelOutput")

  # iterate each output type
  for (modelOutput in modelOutputs) {
    subsetDF <- filter(modelOutputsTable, ModelOutput == modelOutput)
    uniqueFiles <- unique(subsetDF$ImageFile)
    # skip if only NA
    if (length(uniqueFiles)==1 && is.na(uniqueFiles)) {
      message(modelOutput, " not available; skipping.")
      next
    }
  
    # dimensions
    if (length(parameters)>1) {
      ncol <- n_distinct(subsetDF[[parameters[2]]])
      nrow <- n_distinct(subsetDF[[parameters[1]]])
      colNames <- unique(as.character(subsetDF[[parameters[2]]]))
      rowNames <- unique(as.character(subsetDF[[parameters[1]]]))
      colTitle <- parameterCW$CleanName[parameterCW$Name == parameters[2]]
      rowTitle <- parameterCW$CleanName[parameterCW$Name == parameters[1]]
    } else {
      ncol <- n_distinct(subsetDF[[parameters[1]]])
      nrow <- 1
      colNames <- unique(as.character(subsetDF[[parameters[1]]]))
      rowNames <- NULL
      colTitle <- parameterCW$CleanName[parameterCW$Name == parameters[1]]
      rowTitle <- NULL
    }
  
    # read PNG files into raster grobs
    img_grobs <- lapply(subsetDF$ImageFile, function(f) {
      rasterGrob(readPNG(f), interpolate = FALSE)
    })

    # build list: corner + column labels + row labels (if any) + panels
    g_list <- list(nullGrob())
    # column labels
    for (lbl in colNames) {
      g_list <- c(g_list, list(textGrob(lbl, gp = gpar(fontsize = 14))))
    }
    # row labels
    if (!is.null(rowNames)) {
      for (lbl in rowNames) {
        g_list <- c(g_list, list(textGrob(lbl, rot = 90, gp = gpar(fontsize = 14))))
      }
    }
    # panel images
    g_list <- c(g_list, img_grobs)

    # calculate index offsets
    cornerIdx   <- 1
    colLblStart <- cornerIdx + 1
    if (!is.null(rowNames)) {
      rowLblStart <- cornerIdx + ncol + 1
      imgStart <- rowLblStart + nrow
    } else {
      imgStart <- cornerIdx + ncol + 1
      }

    # layout matrix dimensions
    total_rows <- nrow + 1
    total_cols <- ncol + 1
    

    # initialize layout matrix
    layout_mat <- matrix(NA, nrow = total_rows, ncol = total_cols)
    # corner
    layout_mat[1, 1] <- cornerIdx
    # top column labels
    for (c in 2:total_cols) {
      layout_mat[1, c] <- cornerIdx + (c - 1)
    }
    # side row labels
    for (r in 2:total_rows) {
      layout_mat[r, 1] <- cornerIdx + ncol + (r - 1)
    }
    # panels
    for (r in 2:total_rows) {
      for (c in 2:total_cols) {
        layout_mat[r, c] <- imgStart + (r - 2) * ncol + (c - 2)
      }
    }

    # fixed widths/heights for grid
    widths  <- unit.c(unit(0.8, "cm"), rep(unit(1, "null"), ncol))
    heights <- unit.c(unit(0.6, "cm"), rep(unit(1, "null"), nrow))

    # assemble matrix with padding
    matrix_grob <- arrangeGrob(
      grobs = g_list,
      layout_matrix = layout_mat,
      widths = widths,
      heights = heights,
      respect = TRUE,
      padding = unit(4, "pt")
    )

    # overall axis titles
    col_title <- textGrob(
      colTitle,
      gp = gpar(fontsize = 16)
    )
    row_title <- if (!is.null(rowNames)) {
      textGrob(
        rowTitle,
        gp = gpar(fontsize = 16),
        rot = 90
      )
    } else {
      NULL
    }

    # wrap with axis titles
    title_widths  <- unit.c(unit(0.8, "cm"), unit(1, "null"))
    title_heights <- unit.c(unit(0.6, "cm"), unit(1, "null"))

    if (!is.null(rowNames)) {
      top_row <- arrangeGrob(
        nullGrob(), col_title,
        ncol = 2, widths = title_widths
      )
      middle <- arrangeGrob(
        row_title, matrix_grob,
        ncol = 2, widths = title_widths
      )
      final_grob <- arrangeGrob(
        top_row, middle,
        ncol = 1, heights = title_heights
      )
    } else {
      final_grob <- arrangeGrob(
        col_title, matrix_grob,
        ncol = 1, heights = title_heights
      )
    }

    # save to high-res PNG
    outfile <- file.path(outputPath, paste0(modelOutput, "Matrix.png"))
    png(
      filename = outfile,
      width  = (0.8 / 2.54) + ncol * 3,
      height = (0.6 / 2.54) + nrow * 3,
      units   = "in",
      res     = 600
    )
    grid.draw(final_grob)
    dev.off()
  }
}


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