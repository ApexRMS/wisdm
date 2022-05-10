## -------------------------------------
## wisdm - prep spatial multiprocessing
## ApexRMS, May 2022
## -------------------------------------

# built under R version 4.1.1
# this transformer pulls in a template raster and creates a tiling raster for 
# spatial multiprocessing 

# source dependencies ----------------------------------------------------------

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "0-dependencies.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myLibrary <- ssimLibrary()
myProject <- rsyncrosim::project()
myScenario <- scenario()
# datasheet(myScenario)

# determine if transformer needs to be run -------------------------------------

runTransformer <- datasheet(myScenario, "core_Multiprocessing")$EnableMultiprocessing

if(runTransformer){ # if true, run transformer
  
  # path to ssim directories
  ssimTransDir <- ssimEnvironment()$TransferDirectory 
    
  # Read in datasheets
  runControlSheet <- datasheet(myScenario, "RunControl", optional = T)
  templateSheet <- datasheet(myScenario, "TemplateRaster")
  spatialMulitprocessingSheet <- datasheet(myScenario, "corestime_Multiprocessing")
  
  # Set defaults -----------------------------------------------------------------
  
  ## Run control sheet
  if(nrow(runControlSheet)<1){
    runControlSheet <- addRow(runControlSheet, list(1,1,0,0))
    saveDatasheet(myScenario, runControlSheet, "RunControl")
  }

  # Setup multiprocessing ------------------------------------------------------
  
  if(nrow(templateSheet)<1){ stop("ERROR: template raster is missing") }
  
  # load template raster
  templateRaster <- rast(templateSheet$RasterFilePath)
  
  tileData <- writeStart(templateRaster, filename = file.path(ssimTransDir, "temp.tif"), overwrite = T)
  writeStop(templateRaster)
  tileData$row <- c(1, 651) # remove
  tileData$nrows <- c(650, tileData$nrows-650) # remove
  tileData$n <- 2 # remove
  
  if(tileData$n > 1){
    
    tileSize <- nrow(templateRaster)*ncol(templateRaster)/tileData$n
    
    # Calculate dimensions of each tile
    tileHeight <- tileData$nrows[1]
    # tileWidth <- floor(tileSize / tileHeight)
    tileWidth <- ncol(templateRaster) # remove
    
    # Calculate number of rows and columns
    ny <- tileData$n
    # nx <- ceiling(ncol(templateRaster) / tileWidth)
    nx <- 1 # remove
    
    # Generate a string of zeros the width of one tile
    oneTileWidth <- rep(0, tileWidth)
    
    # Generate one line of one row, repeat to the height of one row
    oneRow <-
      as.vector(vapply(seq(nx), function(i) oneTileWidth + i, FUN.VALUE = numeric(tileWidth))) %>%
      `[`(1:ncol(templateRaster)) %>% # Trim the length of one row to fit in template
      rep(tileHeight)
    
    # Write an empty raster with the correct metadata to file
    tileRaster <- rast(templateRaster)
    
    # Write tiling to file row-by-row
    writeStart(tileRaster, filename = file.path(ssimTransDir, "tile.tif"),  overwrite=TRUE)
    for(i in seq(tileData$n)) {
      if(tileData$nrows[i] < tileHeight)
        oneRow <- oneRow[1:(ncol(tileRaster) * tileData$nrows[i])]
        writeValues(x = tileRaster, v = oneRow, start = tileData$row[i], nrows = tileData$nrows[i])
        oneRow <- oneRow + nx
    }
    writeStop(tileRaster)
    
    # save tiling raster to core datasheet ---------------------------------------
    
    spatialMulitprocessingSheet <- addRow(spatialMulitprocessingSheet, 
                                          c(MaskFileName = file.path(ssimTransDir, "tile.tif")))
    
    saveDatasheet(myScenario, spatialMulitprocessingSheet, "corestime_Multiprocessing")
  
    } # end if statement: tileData$n > 1  
} else {
  stop("ERROR: Multiprocessing is not enabled. Enable multiprocessing to prepare a tiling raster.")
} # end run transformer
