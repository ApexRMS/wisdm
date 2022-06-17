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

# Active scenario
myScenario <- scenario()

# path to ssim directories
ssimTransDir <- ssimEnvironment()$TransferDirectory 
    
# Read in datasheets
templateSheet <- datasheet(myScenario, "TemplateRaster")
spatialMulitprocessingSheet <- datasheet(myScenario, "corestime_Multiprocessing")
  

# Setup multiprocessing --------------------------------------------------------
  
  if(nrow(templateSheet)<1){ stop("Template raster is missing") }
  
  # load template raster
  templateRaster <- rast(templateSheet$RasterFilePath)
  tileCount <- templateSheet$TileCount
  
  # if tile count is provided
  if(!is.na(tileCount)){
    
    nRow <- round(nrow(templateRaster)/tileCount, 0)
    startRows <- seq(1,nrow(templateRaster), by = nRow)
    nRows <- rep(nRow, length(startRows))
  
    if(tail(nRows,1)+tail(startRows,1)-1 != nrow(templateRaster)){
      nRows[length(nRows)] <- nrow(templateRaster)-tail(startRows,1)+1
      # sum(nRows) == nrow(templateRaster)
    }
    
    tileData <- list(row = startRows, 
                   nrows = nRows, 
                   n = tileCount)
    
    suggestedTiles <- writeStart(templateRaster, filename = file.path(ssimTransDir, "temp.tif"), overwrite = T)
    writeStop(templateRaster)
    # tempRast <- raster::raster(templateSheet$RasterFilePath)
    # suggestedTiles <- raster::blockSize(tempRast)$n
    # remove(tempRast)
    
    if(abs(suggestedTiles$n-tileData$n)>2){ message(paste0("The output tiling raster includes ", 
                                                      tileData$n, " tiles when the suggested tile count is ",
                                                      suggestedTiles$n, ". This may cause longer processing time."))}
   
    
  } else {
    # To Do: update to terra call once blockSize is supported by terra 
    tileData <- writeStart(templateRaster, filename = file.path(ssimTransDir, "temp.tif"), overwrite = T)
    writeStop(templateRaster)
    
    # tempRast <- raster::raster(templateSheet$RasterFilePath)
    # tileCount <- raster::blockSize(tempRast)$n
    # remove(tempRast)
    # 
    # nRow <- round(nrow(templateRaster)/tileCount, 0)
    # startRows <- seq(1,nrow(templateRaster), by = nRow)
    # nRows <- rep(nRow, length(startRows))
    # 
    # if(tail(nRows,1)+tail(startRows,1)-1 != nrow(templateRaster)){
    #   nRows[length(nRows)] <- nrow(templateRaster)-tail(startRows,1)+1
    #   # sum(nRows) == nrow(templateRaster)
    # }
    # 
    # tileData <- list(row = startRows, 
    #                  nrows = nRows, 
    #                  n = tileCount)
    
    templateSheet$TileCount <- tileData$n
    saveDatasheet( myScenario, templateSheet, "TemplateRaster")
    
  }
  
  if(tileData$n > 1){
    
    tileSize <- round(nrow(templateRaster)*ncol(templateRaster)/tileData$n,0)
    
    # Calculate dimensions of each tile
    tileHeight <- tileData$nrows[1]
    tileWidth <- floor(tileSize / tileHeight)

    # Generate a string of zeros the width of one tile
    oneTileWidth <- rep(0, tileWidth)
    
    # Generate one line of one row, repeat to the height of one row
    oneRow <-
      as.vector(vapply(seq(1), function(i) oneTileWidth + i, FUN.VALUE = numeric(tileWidth))) %>%
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
        oneRow <- oneRow + 1
    }
    writeStop(tileRaster)
    
    # save tiling raster to core datasheet ---------------------------------------
    
    spatialMulitprocessingSheet <- addRow(spatialMulitprocessingSheet, 
                                          c(MaskFileName = file.path(ssimTransDir, "tile.tif")))
    
    saveDatasheet(myScenario, spatialMulitprocessingSheet, "corestime_Multiprocessing")
  
    } else { # if tileData$n == 1  
      message("The default settings suggest a tile count of 1; no tiling raster created. A tile count can be set to overide the default settings.")
    }
  