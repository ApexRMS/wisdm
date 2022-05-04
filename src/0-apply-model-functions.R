## ------------------------------
## wisdm - apply model functions  
## ApexRMS, April 2022       
## ----------------------------- 


# Proc.tiff Function -----------------------------------------------------------
  
  # Description:
  # This function is used to make predictions using a number of .tiff image inputs
  # in cases where memory limitations don't allow the full images to be read in.
  
  proc.tiff <- function(fit.model,
                        modType,
                        mod.vars,
                        raster.files = NULL,
                        output.options = NULL,
                        factor.levels = NA,
                        temp.directory, 
                        multiprocessing.cores = 1,
                        train.dat,
                        pred.fct,
                        tsize=2.0,
                        NAval = -3.399999999999999961272e+38)
    {
    

  # Arguments:
  # vname: names of variables used for prediction.  must be same as filenames and variables
  #   in model used for prediction. do not include a .tif extension as this is added in code.
  # fpath: path to .tif files for predictors. use forward slashes and end with a forward slash ('/').
  # myfun: prediction function.  must generate a vector of predictions using only a
  #   dataframe as input.
  # outfile:  name of output file.  placed in same directory as input .tif files.  should
  #   have .tif extension in name.
  # tsize: size of dataframe used for prediction in MB.  this controls the size of tiles
  #   extracted from the input files, and the memory usage of this function.
  # NAval: this is the NAvalue used in the input files.
  
  # Description:
  # This function reads in a limited number of lines of each image (specified in terms of the
  # size of the temporary predictor dataframe), applies a user-specified
  # prediction function, and stores the results as matrix.  Alternatively, if an
  # output file is specified, a file is written directly to that file in .tif format to
  # the same directory as the input files.  Geographic information from the input images
  # is retained.
  
  # Start of function #
  makeProb <- output.options$MakeProbabilityMap
  makeMESS <- output.options$MakeMessMap 
  makeMOD <- output.options$MakeModMap
 
  if(is.null(factor.levels)){ factor.levels <- NA }
  
  nVars <- length(mod.vars)
  
  if(nVars < 1){
    makeMESS <- FALSE
    message("NOTE: MESS map not generated. Model must have 1 or more predictor variables to make a MESS map.")
  } 
  if(nVars == 1){
    makeMOD <- FALSE 
    message("NOTE: MOD map not generated. Model must have at least 1 predictor variables to make a MOD map.")
  }
  
  # get spatial reference info from existing image file
  RasterInfo <- rast(raster.files[1])
  
  # tilize(templateRaster = raster::raster(raster.files[1]), 
  #        # tempfilename = paste0(ssimTempDir,"tempTile.tif"), 
  #        filename = paste0(ssimTempDir,"tile.tif"),
  #        tileSize = 1e6L)
  
  dims <- dim(RasterInfo)[1:2]
  # ps <- res(RasterInfo)
  # ll <- RasterInfo@ptr$extent$vector[c(1,3)] # lower left origin x,y
  # pref <- crs(RasterInfo)
  
  # setting tile size
  MB.per.row <- dims[2]*nVars*32/8/1000/1024 # TO Do: check why this calculation is used
 
  if(makeMESS){ # use more blocks for mess
    MB.per.row <- MB.per.row*8 
  } 
  nrows <- min(round(tsize/MB.per.row),dims[1])
 
  tr <- list()
  tr$row <- seq(1,dims[1], by = nrows)
  tr$nrows <- diff(tr$row)
  if(tail(tr$row,1)<dims[1]){ tr$nrows <- c(tr$nrows, (dims[1] - tail(tr$row,1)+1)) }
  tr$n <- length(tr$row) 
  
  gi <- rgdal::GDALinfo(raster.files[1])
  if(!is.na(match("AREA_OR_POINT=Point",attr(gi,"mdata")))){                                                                          
    xx<-RasterInfo  # this shifts by a half pixel
    nrow(xx) <- nrow(xx) - 1
    ncol(xx) <- ncol(xx) - 1
    rs <- res(xx)
    xmin(RasterInfo) <- xmin(RasterInfo) - 0.5 * rs[1]
    xmax(RasterInfo) <- xmax(RasterInfo) - 0.5 * rs[1]
    ymin(RasterInfo) <- ymin(RasterInfo) + 0.5 * rs[2]
    ymax(RasterInfo) <- ymax(RasterInfo) + 0.5 * rs[2]
  }
  
  FactorInd <- which(!is.na(match(mod.vars,names(factor.levels))),arr.ind=TRUE)
  if((nVars-length(FactorInd))==0) { # turn this off if only one factor column was selected
    makeMESS <- makeMOD <- FALSE 
  }
    
  # create some temporary folders    
  if(makeProb){ dir.create(file.path(temp.directory,"ProbTiff")) }
  if(makeMESS){ dir.create(file.path(temp.directory,"MESSTiff")) } 
  if(makeMOD){ dir.create(file.path(temp.directory,"ModTiff")) } 
  if(makeResid){ dir.create(file.path(temp.directory,"ResidTiff")) }
  
  # turn off multicore in certain circumstances
  if(is.na(multiprocessing.cores)){ multCore <- FALSE } else { multCore <- TRUE }
  if(tr$n<10 | getRversion()<2.14){ multCore <- FALSE }
  
  if(multCore){
    
    library(parallel)
    
    # identify available cores and assign tiled raster processing 
    tile.start <- seq(from=1, to=tr$n, by=ceiling(tr$n/multiprocessing.cores)) 
    cl <- makeCluster(multiprocessing.cores) 
    parLapply(cl, X = tile.start, fun = parRaster, 
              nToDo = ceiling(tr$n/multiprocessing.cores),
              pred.fct=pred.fct,
              dims=dims,
              tr=tr,
              makeProb=makeProb,
              makeMESS=makeMESS,
              makeMOD=makeMOD,
              RasterInfo=RasterInfo,
              template=RasterInfo,
              nVars=nVars,
              temp.directory = temp.directory,
              raster.files=raster.files,
              mod.vars=mod.vars,
              NAval=NAval,
              factor.levels=factor.levels,
              fit.model=fit.model,
              modType = modType,
              train.dat=train.dat)
    stopCluster(cl)
  } else {  # multicore is slower for small tiffs so we won't do it and the library is not available prior to 2.14
    # also due to multicore multiinstance R issues we're currently only running it on condor or when running synchronously
    parRaster(start.tile=1, 
              nToDo = tr$n,
              dims=dims,
              tr=tr,
              makeProb=makeProb,
              makeMESS=makeMESS,
              makeMOD=makeMOD,
              nVars=nVars,
              temp.directory = temp.directory,
              raster.files=raster.files,
              mod.vars=mod.vars,
              NAval=NAval,
              factor.levels=factor.levels,
              RasterInfo=RasterInfo,
              template=RasterInfo,
              fit.model=fit.model,
              modType = modType,
              pred.fct=pred.fct,
              train.dat=train.dat)
  }
  if(length(l <- list.files(paste0(temp.directory,"ProbTiff\\"),pattern="_prob_map.txt",full.names=TRUE))!= 0){ unlink(l) } 
  return(0)
  }

## Tilize function -------------------------------------------------------------
  
  # Function to generate a tiling mask given template
  # - To avoid holding the entire raster in memory, the raster is written directly
  #   to file row-by-row. This way only one row of the tiling needs to be held in
  #   memory at a time. This is also why the number of rows (and implicitly the size
  #   of a given row) cannot be chosen manually
  # - minProportion is used to determine the size threshold for consolidating tiles
  #   that are too small into neighboring tiles. Represented as a porportion of a
  #   full tile
  # 
  # tilize <- function(templateRaster, 
  #                    filename, 
  #                    tempfilename, 
  #                    tileSize) {
  #   
  #   # Calculate recommended block size of template
  #   blockInfo <- raster::blockSize(templateRaster)
  #   
  #   # Check that the blockSize is meaningful
  #   # - This should only matter for very small rasters, such as in test mode
  #   if(max(blockInfo$nrows) == 1){
  #     blockInfo <- list(row = 1, nrows = nrow(templateRaster), n = 1)
  #   }
  #     
  #   # Calculate dimensions of each tile
  #   tileHeight <- blockInfo$nrows[1]
  #   tileWidth <- floor(tileSize / tileHeight)
  #   
  #   # Calculate number of rows and columns
  #   ny <- blockInfo$n
  #   nx <- ceiling(ncol(templateRaster) / tileWidth)
  #   
  #   # Generate a string of zeros the width of one tile
  #   oneTileWidth <- rep(0, tileWidth)
  #   
  #   # Generate one line of one row, repeat to the height of one row
  #   oneRow <-
  #     as.vector(vapply(seq(nx), function(i) oneTileWidth + i, FUN.VALUE = numeric(tileWidth))) %>%
  #     `[`(1:ncol(templateRaster)) %>% # Trim the length of one row to fit in template
  #     rep(tileHeight)
  #   
  #   # Write an empty raster with the correct metadata to file
  #   tileRaster <- raster::raster(templateRaster)
  #   
  #   # Write tiling to file row-by-row
  #   tileRaster <- raster::writeStart(tileRaster, filename,  overwrite=TRUE)
  #   for(i in seq(blockInfo$n)) {
  #     if(blockInfo$nrows[i] < tileHeight)
  #       oneRow <- oneRow[1:(ncol(tileRaster) * blockInfo$nrows[i])]
  #     #browser()
  #     tileRaster <- writeValues(tileRaster, oneRow, blockInfo$row[i])
  #     oneRow <- oneRow + nx
  #   }
  #   tileRaster <- writeStop(tileRaster)
  #   
  #   # # Mask raster by template
  #   # tileRaster <-
  #   #   maskRaster(inputRaster = tileRaster, filename = tempfilename, maskingRaster = templateRaster)
  #   # 
  #   # # Consolidate small tiles into larger groups
  #   # # - We want a map from the original tile IDs to new consolidated tile IDs
  #   # reclassification <-
  #   #   # Find the number of cells in each tile ID
  #   #   tabulateRaster(tileRaster) %>%
  #   #   # Sort by ID to approximately group tiles by proximity
  #   #   arrange(value) %>%
  #   #   # Consolidate into groups up to size tileSize
  #   #   transmute(
  #   #     from = value,
  #   #     to   = consolidateGroups(freq, tileSize)) %>%
  #   #   as.matrix()
  #   # 
  #   # # Reclassify the tiling raster to the new consolidated IDs
  #   # tileRaster <-
  #   #   reclassify(
  #   #     tileRaster,
  #   #     reclassification,
  #   #     filename = filename,
  #   #     overwrite = T)
  #   
  #   return(tileRaster)
  # } 

### maskRaster function --------------------------------------------------------
  
  # # A memory safe implementation of raster::mask() optimized for large rasters
  # # - Requires an output filename, slower than raster::mask for small rasters
  # # - Input and mask rasters must have same extent (try cropRaster() if not)
  # 
  # maskRaster <- function(inputRaster, filename, maskingRaster){
  #   # Create an empty raster to hold the output
  #   outputRaster <-  raster::raster(inputRaster)
  #   
  #   ## Split output into manageable chunks and fill with data from input
  #   blockInfo <- raster::blockSize(outputRaster)
  #   
  #   ## Calculate mask and write to output block-by-block
  #   outputRaster <- raster::writeStart(outputRaster, filename, overwrite=TRUE)
  #   
  #   # Each block of the mask raster is converted into a multiplicative mask
  #   # and multiplied with the corresponding block of the input raster
  #   for(i in seq(blockInfo$n)) {
  #     maskedBlock <-
  #       raster::getValuesBlock(maskingRaster, row = blockInfo$row[i], nrows = blockInfo$nrows[i]) %>%
  #       maskify %>%
  #       `*`(raster::getValuesBlock(inputRaster, row = blockInfo$row[i], nrows = blockInfo$nrows[i]))
  #     outputRaster <- raster::writeValues(outputRaster, maskedBlock, blockInfo$row[i])
  #   }
  #   
  #   outputRaster <- raster::writeStop(outputRaster)
  #   
  #   return(outputRaster)
  # }
  
#### maskify function ----------------------------------------------------------
  
  # # Function to convert a vector, etc. of number into a multiplicative mask
  # # - Replaces all values with 1, NA remains as NA
  # # - Assumes no negative numbers (specifically, no -1)
  # 
  #  maskify <- function(x) {
  #   x <- x + 1L # Used to avoid divide by zero errors, this is why -1 is not acceptable
  #   return(x / x)
  # }
      
## parRaster function ---------------------------------------------------------
  
  parRaster <- function(start.tile,
                        nToDo,
                        dims,
                        tr,
                        makeProb,
                        makeMESS,
                        makeMOD,
                        nVars,
                        raster.files,
                        temp.directory,
                        mod.vars,
                        NAval,
                        factor.levels,
                        RasterInfo,
                        template, 
                        pred.fct,
                        train.dat,
                        fit.model,
                        modType)
  {
    # packageDir <- Sys.getenv("ssim_package_directory")
    # source(file.path(packageDir, "0-dependencies.R"))
    # source(file.path(packageDir, "0-helper-functions.R"))
    # source(file.path(packageDir, "0-apply-model-functions.R"))
    # library(terra)
    
    # source(paste0(file.path(ScriptPath),"/pred.fct.r",sep=''))
    # source(paste0(file.path(ScriptPath),"/maxent.predict.r",sep=''))
    # source(paste0(file.path(ScriptPath),"/CalcMESS.r",sep=''))
    
    
    # for the last set we have to adjust tr$n based on the number of remaining tiles
    if((start.tile + nToDo) > tr$n){ 
      nToDo <- tr$n - start.tile + 1 
    } 
    
    if(tr$n > (start.tile+nToDo)){
      # hack to change the extent for writing sections to separate files because crop crashes for large files
      start.val <- xyFromCell(RasterInfo, cellFromRowCol(RasterInfo, 
                                                         row = tr$row[(start.tile+nToDo)]-1, 
                                                         col = 1)) - 0.5 * res(RasterInfo)[2]
      end.val <- xyFromCell(RasterInfo, cellFromRowCol(RasterInfo, 
                                                       row = tr$row[start.tile], 
                                                       col = ncol(RasterInfo))) + 0.5 * res(RasterInfo)[2]
      e <- ext(RasterInfo)
      ext(RasterInfo) <- c(e[1:2], start.val[1, 2], end.val[1, 2])
      v <- values(RasterInfo, nrows = as.integer(sum(tr$nrows[start.tile:(start.tile + nToDo - 1)])))
      nrow(RasterInfo) <- as.integer(sum(tr$nrows[start.tile:(start.tile + nToDo - 1)]))
      values(RasterInfo) <- v
    }
    
    outfile.p <- paste0(temp.directory,"ProbTiff\\", modType,"_prob_map.tif")
    outfile.p <- file.path(paste(substr(outfile.p, 1, (nchar(outfile.p) - 4)), ifelse(start.tile == 1, "", start.tile), ".tif", sep = ""))

    # finding a folder that's been created so we can keep the user updated on the progress of creating the files
    if(makeProb){
      outtext <- paste(substr(outfile.p,1,(nchar(outfile.p)-4)),".txt",sep="")
    } else {
      if(makeMESS){
        outtext <- paste0(temp.directory,"MESSTiff\\", modType,"_mess_map.txt")
        if(makeMOD){
          outtext <- paste0(temp.directory,"ModTiff\\", modType,"_mod_map.txt")
        }
      }
    }

    capture.output(cat(paste(nToDo,"tiles to do\n")),file=outtext,append=TRUE)

    for (run in 1:2){
      # start up rasters we need   
      if(makeProb & run == 1){
        continuousRaster <- rast(RasterInfo)
        writeStart(x = continuousRaster, 
                   filename = outfile.p, 
                   overwrite = TRUE,
                   datatype ='INT1U',
                   gdal=c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
      }
      
      if(makeMESS & run == 2){
        MessRaster <- rast(RasterInfo)
        writeStart(MessRaster, 
                   filename = sub("ProbTiff", "MESSTiff", sub("prob", "mess", outfile.p)), 
                   overwrite = TRUE,
                   datatype ='INT2S',
                   gdal=c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
        
        if(makeMOD){
          ModRaster <- RasterInfo
          writeStart(ModRaster, 
                     filename = sub("ProbTiff", "ModTiff", sub("prob", "MoD", outfile.p)), 
                     overwrite = TRUE,
                     datatype ='INT1U', 
                     gdal=c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
        }
        
        # order the training data so that we can consider the first and last row  only in mess calculations
        train.dat <- train.dat[ ,match(mod.vars, names(train.dat))]
        for(k in 1:nVars){ train.dat[ ,k] <- sort(train.dat[ ,k]) } 
      }
      
        HasTemplate <- FALSE
        TemplateMask <- NA
        templateRast <- template
      
      # if(class(templateRast)=="try-error"){ #so that we can move a session folder to a new computer
      #   template <- file.path(dirname(maDir), basename(template))
      #   templateRast <- try(raster(template), silent = TRUE) 
      #   if(class(templateRast) == "try-error") HasTemplate = FALSE
      # }
      
      for (i in start.tile:min(start.tile + nToDo - 1, length(tr$row))){
        # i <- 1
        capture.output(cat(paste("starting tile", i, Sys.time(), "\n")), file = outtext, append = TRUE)
        
        # alter the write start location because we always start at position 1                                   
        writeLoc <- ifelse((start.tile - 1) == 0, tr$row[i], tr$row[i] - sum(tr$nrows[1:(start.tile-1)]))
        
        if(HasTemplate){
          TemplateMask <- values(templateRast, row = tr$row[i], nrows = tr$nrows[i])
  
          if(all(is.na(TemplateMask))){
  
            # if the template is completely NA values, don't read in any other data
            temp <- rep(NA, times = tr$nrow[i] * dims[2])
  
            if(makeMESS & run == 2){ pred.rng <- rep(NA,length(temp)) }
          }
        }
  
        if(!HasTemplate | !all(is.na(TemplateMask))){
          
          temp <- data.frame(matrix(ncol = nVars, nrow = tr$nrows[i] * dims[2]))
          
          # Setting the first predictor equal to NA where ever the mask is NA
          # fill temp data frame
          for(k in 1:nVars){
            rast_k <- rast(raster.files[k])
            temp[,k] <- as.vector(values(rast_k, row = tr$row[i], nrows = tr$nrows[i]))
          }
          # so we won't write out predictions here
          
          names(temp) <- mod.vars
          
          if(HasTemplate){ temp[is.na(TemplateMask),] <- NA }
          
          # replace missing values 
          temp[temp == NAval] <- NA
          for(m in mod.vars){ temp[is.nan(temp[,m]),m] <- NA }
          
          if(makeMESS & run == 2){
            pred.rng <- rep(NA, nrow(temp))
            names(pred.rng) <- NA
            
            if(any(complete.cases(temp))){
              MessVals <- CalcMESS(rast = temp[complete.cases(temp), ], train.dat = train.dat)
              pred.rng[complete.cases(temp)] <- MessVals[ ,2]
              names(pred.rng)[complete.cases(temp)] <- MessVals[ ,1]
            }
            
            writeValues(x = MessRaster, v = round(pred.rng*100,0), start = writeLoc, nrows = tr$nrows[i])
            if(is.null(names(pred.rng))) names(pred.rng) <- NA
            if(makeMOD) { writeValues(x = ModRaster, v = names(pred.rng), writeLoc, nrows = tr$nrows[i]) }
          }
          
          if(length(mod.vars) == 1){ names(temp) <- mod.vars }
           
          if(!is.na(factor.levels)){
            factor.cols <- match(names(factor.levels), names(temp))
            if(sum(!is.na(factor.cols)) > 0){
              for(j in 1:length(factor.cols)){
                if(!is.na(factor.cols[j])){
                  temp[,factor.cols[j]] <- factor(temp[, factor.cols[j]], levels = factor.levels[[j]]$number, labels = factor.levels[[j]]$class)
                }
              }
            }
          } 
        } 
        
        # will not calculate predictions if all predictors in the region are na
        if(makeProb & run == 1){ 
        
          ifelse(sum(complete.cases(temp)) == 0,  
                 preds <- matrix(data = NA, nrow = dims[2], ncol = tr$nrows[i]),
                 preds <- t(matrix(pred.fct(model = fit.model, x = temp), ncol = dims[2], byrow = T)))
          
          preds <- round((preds*100), 0)
          
          writeValues(x = continuousRaster, v = preds, start = writeLoc, nrows = tr$nrows[i]) 
        }
      } #end of the big for loop
      
      # end.seq <- c(tr$row, dims[1] + 1)
      
      if(makeProb & run == 1){ writeStop(continuousRaster) }
      if(makeMESS & run == 2){ writeStop(MessRaster) }
      if(makeMOD){ 
        writeStop(ModRaster)
        d <- data.frame(as.integer(seq(1:length(train.dat))), names(train.dat))
        names(d) = c("Value","Class")
        # ModRaster@file@datanotation <- "INT1U"
        write.dbf(d, sub(".tif", ".tif.vat.dbf", ModRaster@file@name), factor2char = TRUE, max_nchar = 254)
      }
    } # end run loop
  }

### Calculate MESS Map function ------------------------------------------------
  
  CalcMESS <- function(rast,train.dat){  
    
    # if anything is out of range, return it before calculating sums                       
    min.train <- train.dat[1,] #because we sorted
    max.train <- train.dat[nrow(train.dat),]
    output <- data.frame(matrix(data=NA,nrow=nrow(rast),ncol=ncol(rast)))
    for(k in 1:length(min.train)){
      output[,k] <- my.min(rast.val=rast[,k], as.numeric(min.train[,k]), as.numeric(max.train[,k]))
    }
    Indx <- apply(output,1,which.min)
    output <- apply(output,1,min)
    if(sum(output>0)==0){ return(data.frame(Indx=Indx,output=output)) }
    
    f <- data.frame(matrix(dat=NA,nrow=sum(output>0),ncol=ncol(rast)))
    
    for(k in 1:length(min.train))  {
      f[,k]<-100*mapply(vecSum,rast[output>0,k],MoreArgs=list(vect=train.dat[,k]))/nrow(train.dat)
      f[,k]<-I(f[,k]<=50)*2*f[,k]+I(f[,k]>50)*2*(100-f[,k])
    }
    
    Indx[output>0] <- apply(f,1,which.min)
    f <- apply(f,1,min)
    output[output>0] <- f
    return(data.frame(Indx=Indx,output=output))
  }
  
#### My Min Function -----------------------------------------------------------
  
  my.min <- function(rast.val,min.train,max.train){
    pmin((rast.val-min.train),(max.train-rast.val))/(max.train-min.train)*100
  }

#### Vector Sum function -------------------------------------------------------

  vecSum <- function(v,vect){ sum(v > vect) }
    
### Predict Surface function ---------------------------------------------------
  
  Pred.Surface <- function(object, 
                           model, 
                           filename="", 
                           na.rm=TRUE,
                           NAval) {
    
    predrast <- rast(object)
    # filename <- trim(filename)
    firstrow <- 1
    firstcol <- 1
    ncols <- ncol(predrast)
    lyrnames <- names(object)
    xylyrnames <- c('x', 'y', lyrnames)
    v <- matrix(NA, ncol=nrow(predrast), nrow=ncol(predrast))
    na.rm <- FALSE
    
    tr <- blockSize(predrast, n=nlayers(object)+5)
    ablock <- 1:(ncol(object) * tr$nrows[1])
    napred <- rep(NA, ncol(predrast)*tr$nrows[1])
    predrast <- writeStart(predrast, filename=filename,overwrite=TRUE)
    ############################################################
    for (i in 1:tr$n) {
      if (i==tr$n) { 
        ablock <- 1:(ncol(object) * tr$nrows[i])
        napred <- rep(NA, ncol(predrast) * tr$nrows[i])
      }
      rr <- firstrow + tr$row[i] - 1
      p <- xyFromCell(predrast, ablock + (tr$row[i]-1) * ncol(predrast)) 
      p <- na.omit(p)
      blockvals <- data.frame(x=p[,1], y=p[,2])
      if (na.rm) {
        blockvals <- na.omit(blockvals)		
      }
      if (nrow(blockvals) == 0 ) {
        predv <- napred
      } else {
        
        predv <- predict(model, blockvals)
        predv[is.na(predv)]<-NA
      }
      if (na.rm) {  
        naind <- as.vector(attr(blockvals, "na.action"))
        if (!is.null(naind)) {
          p <- napred
          p[-naind] <- predv
          predv <- p
          rm(p)
        }
      }
      
      # to change factor to numeric; should keep track of this to return a factor type RasterLayer
      predv = as.numeric(predv)
      predrast <- writeValues(predrast, predv, tr$row[i])
      #NAvalue(predrast)<-NAval
      print(i)
    }
    
    predrast <- writeStop(predrast)
    
    return(predrast)
  }
  