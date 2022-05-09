## -------------------------------------
## wisdm - prep spatial multiprocessing
## ApexRMS, May 2022
## -------------------------------------

# built under R version 4.1.1
# this transformer pulls in a template raster and creates a tiling raster for 
# spatial mulitprocessing 

# source dependencies ----------------------------------------------------------

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "0-dependencies.R"))
# source(file.path(packageDir, "0-helper-functions.R"))
# source(file.path(packageDir, "0-apply-model-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myLibrary <- ssimLibrary()
myProject <- rsyncrosim::project()
myScenario <- scenario()
# datasheet(myScenario)

# path to ssim directories
ssimTransDir <- ssimEnvironment()$TransferDirectory 
  
# Read in datasheets
runControlSheet <- datasheet(myScenario, "RunControl", optional = T)
multiprocessingSheet <- datasheet(myScenario, "core_Multiprocessing")
templateSheet <- datasheet(myScenario, "TemplateRaster")
spatialMulitprocessingSheet <- datasheet(myScenario, "corestime_Multiprocessing")

# Set defaults -----------------------------------------------------------------

## Run control sheet
if(nrow(runControlSheet)<1){
  runControlSheet <- addRow(runControlSheet, list(1,1,0,0,T))
  saveDatasheet(myScenario, runControlSheet, "RunControl")
}

# determine if transformer needs to be run -------------------------------------
runTransformer <- multiprocessingSheet$EnableMultiprocessing

if(runTransformer){ # if true, run transformer

  # Setup multiprocessing ------------------------------------------------------

  # load template raster
  templateRaster <- rast(templateSheet$RasterFilePath)
  
  tileData <- writeStart(templateRaster, filename = file.path(ssimTransDir, "temp.tif"), overwrite = T)
  writeStop(templateRaster)
  
  if(tileData$n > 1){
    
    tileSize <- nrow(templateRaster)*ncol(templateRaster)/tileData$n
    
    # Calculate dimensions of each tile
    tileHeight <- tileData$nrows[1]
    tileWidth <- floor(tileSize / tileHeight)
    
    # Calculate number of rows and columns
    ny <- tileData$n
    nx <- ceiling(ncol(templateRaster) / tileWidth)
    
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
} # end run transformer
