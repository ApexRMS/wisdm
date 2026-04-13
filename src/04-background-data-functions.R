## -------------------------
## wisdm - background (pseudo-absence) data functions
## ApexRMS, July 2025
## -------------------------

# Background surface generation function ---------------------------------------
backgroundSurfacePointGeneration <- function(
    sp,         # species
    n,          # number of pseudoabsence points to generate
    template,   # template raster
    outputDir,  # output directory
    dat,        # field data 
    method)     # includes method:('kde', 'mcp'); surface: ('continuous', 'binary'); isopleth
{ # start function  
  # TWO methods:
  #     KDE or MCP
  # TWO options for KDE:
  # continuous or binary (random within isopleth)
  # ONE option for MCP:
  # binary (random within isopleth)
  # if KDE & continous: generate pts based on KDE surface only
  # if KDE & binary: generate points randomly within isopleth given
  # if MCP: generate points randomly within MCP, MCP is generaged by isopleth given
  
  ### ### ### ###
  ### objects that all methods need
  ### ### ### ###
  # terra::terraOptions(tempdir=outputDir)
  
  
  ### ### ### ###
  ### objects that all methods need
  ### ### ### ###
  ### sf vector object of observations. will need to set CRS if obs are different
  occ_sf <- sf::st_as_sf(
    dat[dat$Response==1, c("X", "Y")],
    coords = c("X", "Y"),
    crs = crs(template),
    remove = FALSE
  )
  
  if ('kde' %in% method) { 
    ### set outnames
    kde_bg_out = gsub('/', '\\\\', paste0(outputDir, '/', sp, '_kde_bg_surface.tif'))
    kde_pts_out = paste0(outputDir, '/', sp, '_kde_bg_pts.csv')
    bg_pts_out = paste0(outputDir, '/', sp, '_kde_bg_pts.shp')
    ### don't waste time generating points if they exist
    if(file.exists(kde_pts_out)){
      stop("KDE background points already exist!\nDelete current points to proceed.")
    }
    
    ### get cell resolution of template
    template_res = terra::res(template)[1]
    ### get bounding box of observations and buffer by 10 grid cells
    bbox <- sf::st_bbox(c(
      xmin = min(dat[,"X"]),
      ymin = min(dat[,"Y"]),
      xmax = max(dat[,"X"]),
      ymax = max(dat[,"Y"])),
      crs = crs(template))
    ### get cell width to use for kernel density. cell width determined here is approximately what's needed to create:
    # 2,500 (50 x 50) 
    # 10,000 (100 x 100) 
    # 1,000,000 (1000 x 1000) 
    # tiles in bounding box of observations for generating the kernel density estimate. This keeps resolution reasonable, but still fast. 
    cell_width <- mean(c(
      ((bbox$xmax - bbox$xmin) / 50),
      ((bbox$ymax - bbox$ymin) / 50)
    ))
    ### expand bounding box by 1/2 width and length
    buffer_width <- mean(c(
      ((bbox$xmax - bbox$xmin) / 2),
      ((bbox$ymax - bbox$ymin) / 2)
    ))
    bbox <- bbox |>
      sf::st_as_sfc() |>
      sf::st_buffer(dist = buffer_width)|> 
      sf::st_bbox()
    ### create owin
    win <- spatstat.geom::as.owin(bbox)
    ### set window
    OBSppp <- spatstat.geom::ppp(dat[,"X"], dat[,"Y"], window=win)
    # plot(OBSppp)
    ### set bandwidth for KDE. This should be a tunable parameter. sigma, estimated here, is ~ h/2, where h is the parameter estimated by kernealHD with the href method. 
    bw <- spatstat.explore::bw.ppl(OBSppp) |> round() * 2 # sigma = h/2, where h is estimated from adehabitatHR:kernelUD (original code)
    npts <- spatstat.geom::npoints(OBSppp)
    ### generate KDE layer
    density_ppp <- spatstat.explore::density.ppp(
      OBSppp,
      sigma=bw,
      weights = rep(1 / npts, npts), # ensure that intensity integrates to 1
      eps=cell_width)
    r <- terra::rast(density_ppp)
    terra::crs(r) <- terra::crs(template)
    
    ### rescale weights. this is time consuming for large rasters.
    terra::setMinMax(r, force=TRUE); mm=terra::minmax(r)
    r <- (r-mm[1])/(mm[2]-mm[1]) * 100
    r <- round(r)
    # set all values that are 0 to be NA 
    r <- ifel(r==0, NA, r)
    
    ### create new template that is cropped to sampling grid
    template_small <- terra::crop(template, r)
    ### resample to template. this is time consuming for large rasters.
    r <- terra::resample(r, template_small, method="bilinear")
    ### crop & mask to template
    r <- terra::mask(r, template_small)
    
    # ### write out raster
    # # global(!is.na(r), "sum", na.rm=TRUE) 
    # terra::writeRaster(x = r,
    #                    filename = kde_bg_out,
    #                    overwrite=TRUE,
    #                    gdal=c('COMPRESS=LZW', 'BIGTIFF=YES', 'TILED=YES')) 
    ### get cells that are observations
    occ_cellIDs <- terra::extract(r, vect(occ_sf), cells = TRUE)
    
    ### IF CONTINOUS, SAMPLE WITH WEIGHTS
    if(tolower(method$surface) == "continuous"){
      pts <- data.table::data.table()
      system.time({
        while(nrow(pts) < n){
          ### take spatially random sample all cells in OCC region
          rnd_sample <- terra::spatSample(
            x=r,
            size=n*2,
            method='random',
            cells=TRUE,
            xy=FALSE,
            replace=TRUE,
            as.raster=FALSE,
            as.points=FALSE,
            na.rm=TRUE) |> data.table::data.table()
          names(rnd_sample)[1] = "cell"
          names(rnd_sample)[2] = "weights"
          
          ### drop cells that are occurrence
          rnd_sample[rnd_sample$cell%in%occ_cellIDs$cell, "weights"] <- 0
          
          ### sample the sample with weights 
          selected_cells <- sample(
            rnd_sample$cell,
            size = 1000,
            replace = TRUE,
            prob = rnd_sample$weights
          )
          
          ### combine sample set with complete set, drop duplicates if any
          pts <- rbind(pts, data.table::as.data.table(xyFromCell(r, selected_cells)))
          pts <- unique(pts)
          pts <- pts[complete.cases(pts), ]
        }
      }) 
      ### reduce pts to size of background points requested. random sample is all that's needed here. not sampling by prob. Just thinning to limit size
      if(nrow(pts) > n){
        pts_ind <- sample(nrow(pts), size=n, replace=FALSE)
        pts <- pts[pts_ind,]
      }
      ### write out
      names(pts) <- c("X", "Y")
      pts$Response <- -9998
      data.table::fwrite(pts, kde_pts_out)
      # occ_sf <- sf::st_as_sf(
      #   pts[, c("X", "Y")],
      #   coords = c("X", "Y"),
      #   crs = crs(template),
      #   remove = FALSE
      # ) |> terra::vect()
      # terra::writeVector(occ_sf, bg_pts_out, overwrite = T)
      
    } # end tolower(method$surface) == "continuous"
    
    ### if BINARY, SAMPLE WITH OUT WEIGHTS
    if(tolower(method$surface) == 'binary'){
      ### set outnames
      # Normalize to create a probability surface
      total_sum <- terra::global(r, "sum", na.rm = TRUE)[1,1]
      f <- r / total_sum
      # Get values and cell numbers for non-NA cells only
      v <- terra::values(f, mat = FALSE)
      non_na_cells <- which(!is.na(v))
      vals <- v[non_na_cells]; rm(v); gc()
      # Get coordinates of non-NA cells using cell numbers
      coords <- terra::xyFromCell(f, non_na_cells); rm(non_na_cells); gc()
      # Sort values and coords by decreasing probability
      ord <- order(vals, decreasing = TRUE)
      vals_sorted <- vals[ord]; rm(vals); gc()
      # coords_sorted <- coords[ord, ]; rm(coords); gc()
      # Cumulative sum to find 95% threshold
      cum_prob <- cumsum(vals_sorted)
      isopleth <- method$isopleth/100
      idx <- which(cum_prob >= isopleth)[1]
      threshold <- vals_sorted[idx]; rm(idx, vals_sorted); gc()
      # Create binary mask of cells >= threshold
      outRast <- terra::ifel(r >= threshold, 1, NA)
      ### mask to template to get holes
      outRast <- mask(outRast, template_small)
      # ### write out raster
      # terra::writeRaster(x <- outRast,
      #                    filename = kde_bg_out,
      #                    overwrite=TRUE,
      #                    gdal=c('COMPRESS=LZW', 'BIGTIFF=YES', 'TILED=YES'))
      
      # plot(template)
      # plot(outVect, add=TRUE)
      # plot(pts, add=TRUE)
      
      ### get CELL IDs of observations
      occ_cellIDs <- terra::extract(outRast, vect(occ_sf), cells = TRUE)
      
      ### RANDOMLY SAMPLE RASTER WITH NO WEIGHTS
      ### then remove occ locations & dups and sample again. 
      pts <- data.table::data.table()
      while(nrow(pts) < n){
        ### take spatially random sample of entire study area
        system.time({
          rnd_sample <- terra::spatSample(
            outRast,
            size=n*2,
            method='random',
            cells=TRUE,
            xy=TRUE,
            replace=FALSE,
            as.raster=FALSE,
            as.points=FALSE,
            na.rm=TRUE) |> data.table::data.table()
          names(rnd_sample)[1] = "cell"
          names(rnd_sample)[2] = "X"
          names(rnd_sample)[3] = "Y"
        }) # 
        
        ### drop cells that are occurrence
        rnd_sample <- rnd_sample[!rnd_sample$cell%in%occ_cellIDs$cell, ] 
        
        ### combine sample set with complete set, drop duplicates if any
        pts <- rbind(pts, rnd_sample)
        pts <- unique(pts)
        pts <- pts[complete.cases(pts), ]
      }
      ### reduce pts to size of background points requested. random sample is all that's needed here. not sampling by prob. Just thinning to limit size
      if(nrow(pts) > n){
        pts_ind <- sample(nrow(pts), size=n, replace=FALSE)
        pts <- pts[pts_ind,]
      }
      ### write out
      pts$Response <- -9998
      data.table::fwrite(pts, kde_pts_out)
      # occ_sf <- sf::st_as_sf(
      #   pts[, c("X", "Y")],
      #   coords = c("X", "Y"),
      #   crs = crs(template),
      #   remove = FALSE
      # ) |> terra::vect()
      # terra::writeVector(occ_sf, bg_pts_out, overwrite = T)
    } # end tolower(method$surface) == 'binary'
  } # end 'kde' %in% method
  
  if("mcp" %in% method){
    ### set outnames
    mcp_pts_out <- paste0(outputDir, '/', sp, '_mcp_bg_pts.csv')
    bg_pts_out <- paste0(outputDir, '/', sp, '_mcp_bg_pts.shp')
    ### don't waste time generating points if they exist
    if(file.exists(mcp_pts_out)){
      stop("MCP background points already exist!\nDelete current points to proceed.")
    }
    
    # get X,Y coords 
    xy <- fieldDataSheet[,c("X", "Y")]
    n_vals <- nrow(xy)
    if (n_vals < 3L){
      stop("Need at least 3 points for a polygon.")
    }
    
    # center points to calc isopleth
    cx <- mean(xy[, 1])
    cy <- mean(xy[, 2])
    
    # Squared distances (no sqrt: avoids extra cost)
    d2 <- (xy[, 1] - cx)^2 + (xy[, 2] - cy)^2
    rm(cx, cy)
    
    # threshold with quantile from isopleth 
    thr <- stats::quantile(d2, probs = method$isopleth / 100, names = FALSE, type = 7)
    keep <- d2 <= thr
    rm(d2, thr)
    
    if (sum(keep) < 3L) stop("Fewer than 3 points remain after trimming.")
    
    # Convex hull indices on trimmed set
    id_trim <- which(keep)
    hull_local <- grDevices::chull(xy[keep, 1], xy[keep, 2])
    hull_idx <- id_trim[hull_local]
    
    # Close polygon ring
    poly_coords <- rbind(xy[hull_idx, , drop = FALSE], xy[hull_idx[1], , drop = FALSE])
    
    # convert to polygon
    polygons <- poly_coords |>
      sf::st_as_sf(coords = c("X", "Y"), crs = crs(template)) |> 
      summarize(geometry = sf::st_combine(geometry)) |> 
      sf::st_cast("POLYGON")
    
    ### set path for bg surface
    bg_out <- paste0(outputDir, '/', sp, '_mcp_bg_surface.tif')
    
    ### convert to vect & set crs
    outVect <- terra::vect(polygons); rm(polygons)
    crs(outVect) <- terra::crs(template)
    
    ### create new template that is cropped to sampling grid
    template_small <- terra::crop(template, outVect)
    
    ### convert to raster
    outRast <- rasterize(outVect, template_small)
    
    ### mask to template to get holes
    outRast <- mask(outRast, template_small)
    
    # ### write out raster
    # terra::writeRaster(x <- outRast,
    #                    filename = bg_out,
    #                    overwrite=TRUE,
    #                    gdal=c('COMPRESS=LZW', 'BIGTIFF=YES', 'TILED=YES'))
    
    # plot(template)
    # plot(outVect, add=TRUE)
    # plot(pts, add=TRUE)
    
    ### get CELL IDs of observations
    occ_cellIDs <- terra::extract(outRast, vect(occ_sf), cells = TRUE)
    
    ### RANDOMLY SAMPLE RASTER WITH NO WEIGHTS
    ### then remove occ locations & dups and sample again. 
    pts <- data.table::data.table()
    while(nrow(pts) < n){
      ### take spatially random sample of entire study area
      system.time({
        rnd_sample <- terra::spatSample(
          outRast,
          size=n*2,
          method='random',
          cells=TRUE,
          xy=TRUE,
          replace=FALSE,
          as.raster=FALSE,
          as.points=FALSE,
          na.rm=TRUE) |> data.table::data.table()
        names(rnd_sample)[1] = "cell"
        names(rnd_sample)[2] = "X"
        names(rnd_sample)[3] = "Y"
      }) # 
      
      ### drop out cells that are occurrence
      rnd_sample <- rnd_sample[!rnd_sample$cell%in%occ_cellIDs$cell, ] 
      
      ### combine sample set with complete set, drop duplicates if any
      pts <- rbind(pts, rnd_sample)
      pts <- unique(pts)
      pts <- pts[complete.cases(pts), ]
    }
    ### reduce pts to size of background points requested. random sample is all that's needed here. not sampling by prob. Just thinning to limit size
    if(nrow(pts) > n){
      pts_ind <- sample(nrow(pts), size=n, replace=FALSE)
      pts <- pts[pts_ind,]
    }
    ### write out
    pts$Response <- -9998
    data.table::fwrite(pts, mcp_pts_out)
    # occ_sf <- sf::st_as_sf(
    #   pts[, c("X", "Y")],
    #   coords = c("X", "Y"),
    #   crs = crs(template),
    #   remove = FALSE
    # ) |> terra::vect()
    # terra::writeVector(occ_sf, bg_pts_out, overwrite = T)
  } # end "mcp" %in% method
  
} # end backgroundSurfacePointGeneration

