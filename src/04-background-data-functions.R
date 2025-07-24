## -------------------------
## wisdm - background (pseudo-absence) data functions
## ApexRMS, July 2025
## -------------------------

# Background surface generation function ---------------------------------------

# library(sp)
# library(adehabitatHR)
# library(spatstat.geom)
# library(sf)

backgroundSurfaceGeneration <- function(sp,         # species
                                        template,   # template raster
                                        mask,       # spatVector object bounding the study area
                                        outputDir,  # output directory
                                        dat,        # field data 
                                        method)     # c('kde', 'mcp')
  { # start function  
    
    xy  <- sp::SpatialPoints(dat[, c('X', 'Y')])
    ud  <- adehabitatHR::kernelUD(xy, extent = 0.5, grid = 150)
    mm  <- ud@coords[order(ud@coords[, 1]), ]
    mm  <- order(ud@coords[, 2])
    m   <- ud@data$ud[mm]
    
    if ('kde' %in% method) { 
      if(tolower(method$surface) == "continuous"){
        kde_bg_out <- gsub('/', '\\\\', paste0(outputDir, '/', sp, '_kde_bg_surface.tif'))
        kde.mat <- matrix(m, nrow = length(unique(ud@coords[, 1])))
      
        t <- terra::rast(ud)
        ext(t) <- ext(ud)
        crs(t) <- crs(template)
        t <- t / max(kde.mat) * 100 # normalize from 0-100
        t <- terra::crop(t, ext(template))
        t <- terra::resample(x = t, y = template, method = 'near')
        t <- terra::ifel(!is.na(template), t, NA)
        
        writeRaster(x = t,
                    filename = kde_bg_out,
                    gdal = c('BIGTIFF=YES', 'TILED=YES', 'COMPRESS=LZW'),
                    overwrite = T )
      }
      if(tolower(method$surface) == 'binary'){
        ver <- adehabitatHR::getverticeshr(ud, method$isopleth)  
        bg_out <- paste0(outputDir, '/', sp, '_kde_bg_surface.shp')
      }
    }
    
    if("mcp" %in% method){
      ver <- adehabitatHR::mcp(xy, percent = method$isopleth)
      bg_out <- paste0(outputDir, '/', sp, '_mcp_bg_surface.shp')
    }
    
    if ("mcp" %in% method | tolower(method$surface) == 'binary') { 
      
      WindowList <- list()
      for (j in 1:length(ver@polygons[[1]]@Polygons)) {
        x <- y <- vector()
        x <- ver@polygons[[1]]@Polygons[[j]]@coords[, 1]
        y <- ver@polygons[[1]]@Polygons[[j]]@coords[, 2]
        xy <- paste(x, y, sep = "")
        x <- x[!(duplicated(xy))]
        y <- y[!(duplicated(xy))]
        
        MyWindow <-
          try(spatstat.geom::owin(poly = list(x = x[length(x):1], y = y[length(x):1])), silent =
                TRUE)
        if (class(MyWindow) == "try-error")
          MyWindow <-
          spatstat.geom::owin(poly = list(x = x[1:length(x)], y = y[1:length(x)]))
        ifelse(j == 1, WindowList <- MyWindow, {
          ifelse(
            spatstat.geom::is.subset.owin(MyWindow, WindowList),
            WindowList <- spatstat.geom::setminus.owin(WindowList, MyWindow),
            WindowList <- spatstat.geom::union.owin(MyWindow, WindowList)
          )
        })
      }
      
      outVect <- terra::vect(sf::st_as_sf(WindowList))
      crs(outVect) <- terra::crs(template)
      outVect <- terra::intersect(outVect, mask)
      terra::writeVector(outVect, bg_out, overwrite = T)
    }
    
  } # end function


# Background point (psuedo-absence) generation function -------------------------

# library(data.table)


backgroundPointGeneration <- function(sp,                           # species name
                                      n,                            # number of pseudoabsence points to generate
                                      method,                       # c('kde','mcp')
                                      outputDir,                    # output directory 
                                      target_file,                  # path to csv with growth form-specific points (target background point pool)
                                      overwrite)                    # overwrite existing files
  { # start function 
  
  if('kde' %in% method){
    
    kde_pts_out <- paste0(outputDir, '/', sp, '_kde_bg_pts.csv')
    
    if(file.exists(kde_pts_out)){
      stop("KDE background points already exist!\nDelete current points to proceed.")
    }
    
    if(tolower(method$surface) == "continuous"){
      ### Find KDE raster
      r <- list.files(path = outputDir, pattern = 'kde.*.tif$', full.names = T)
      r <- terra::rast(r)
      
      valid <- data.table::data.table()
      
      # you can extract many values, very quickly with terra
      inc <- 10000000
      total <- terra::ncell(r)
      
      # since we're not ever going to put a point on the kde surface where values
      # are < 1, surfaces can be subset to the *truly* available values. For 
      # species with very limited distributions, this dramatically decreases the
      # amount of time the code spends searching for valid bg locations.
      
      if(total > inc){
        
        cat('Finding viable cells...\n')
        
        pb <- txtProgressBar(0, 100, style = 3)
        
        for(i in seq(1,ncell(r),inc)){
          
          tmp <- terra::extract(r, i:(i + inc))
          tmp$id <- c(i:(i + inc))
          tmp$kde_bg_surface <- ifelse(tmp[ ,1] < 1, NA, tmp[ ,1])
          tmp <- na.omit(tmp)
          
          valid <- rbind(valid, tmp)
          
          perc_complete <- ((i + inc) / total) * 100
          
          setTxtProgressBar(pb, ifelse(perc_complete <= 100, perc_complete, 100))
          
        }
        
        close(pb)
        
      } else {
        
        tmp <- terra::extract(r, 1:inc)
        tmp$id <- 1:inc
        tmp$kde_bg_surface <- ifelse(tmp[ ,1] < 1, NA, tmp[ ,1])
        tmp <- na.omit(tmp)
        
        valid <- rbind(valid, tmp)
        
      }
      
      pts <- data.table::data.table()
      used <- data.table::data.table()
      
      while(nrow(pts) < n){
        
        tmp <- valid[sample(nrow(valid), 1000),]
        tmp <- tmp[!(tmp$id %in% used$id), ]
        
        if(nrow(tmp) == 0) next()
        
        used <- rbind(used, tmp)
        
        # generate column of random integers between 1 and 100. 
        # If the id field is >= this number, it is retained.
        tmp[, test := as.numeric(sample(1:100, nrow(tmp), replace = T))]
        
        # keep true values only
        tmp <- tmp[, retain := ifelse(tmp[,1] >= test, T, F)][retain == T]
        
        if(nrow(tmp) > 0){
          pts <- rbind(pts, tmp)
        } else {
          next()
        }
      }
      
      if(nrow(pts) > n){
        pts <- pts[sample(1:nrow(pts), n, prob = pts[[1]], replace = F), ]
      }
      
      # convert to spatVector and write out
      v <- as.data.frame(terra::xyFromCell(r, pts[[2]]))
      v$Response <- c(rep(-9998, n))
      colnames(v) <- c("X", "Y", "Response")
      write.csv(v, kde_pts_out, row.names = F)
      
      cat(nrow(v), 'pseudoabsence points generated for continuous KDE surface\n')
    }
    
    if(tolower(method$surface) == "binary"){
      # v <- sf::st_read(list.files(path = outputDir, pattern = 'kde.*shp', full.names = T), quiet = TRUE)
      #
      # if(missing(target_file)) stop('Missing path to target bg dataset!')
      # pts <- readr::read_csv(target_file, show_col_types = F, progress = F) |>
      #   dplyr::select(1:2) |>
      #   sf::st_as_sf(coords = c('X','Y'), crs = sf::st_crs(v))
      # 
      # test <- sf::st_intersects(v, pts, sparse = T)
      # pts <- pts[unlist(test),] |> unique()
      # 
      # if(nrow(pts) > n){
      #   pts <- pts[sample(nrow(pts), n), ]
      # }
      # 
      # pts <- as.data.frame(sf::st_coordinates(pts))
      
      v <- vect(list.files(path = outputDir, pattern = 'kde.*shp', full.names = T))
      pts <- crds(spatSample(v, n, "random"))
      
      pts <- as.data.frame(pts)
      names(pts) <- c("X", "Y")
      pts$Response <- -9998
      write.csv(pts, kde_pts_out, row.names = F)
      
      cat(nrow(pts), 'pseudoabsence points generated for binary KDE surface\n')
    }
  }
  
  if('mcp' %in% method){ 
    
    mcp_pts_out <- paste0(outputDir, '/', sp, '_mcp_bg_pts.csv')
    
    v <- vect(list.files(path = outputDir, pattern = 'mcp.*shp', full.names = T))
    pts <- crds(spatSample(v, n, "random"))
    
    pts <- as.data.frame(pts)
    names(pts) <- c("X", "Y")
    pts$Response <- -9998
    write.csv(pts, mcp_pts_out, row.names = F)
    
    cat(nrow(pts), 'pseudoabsence points generated for MCP surface\n')
  }
  
} # End Function