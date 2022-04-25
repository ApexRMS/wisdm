## ------------------------------
## wisdm - apply model functions  
## ApexRMS, April 2022       
## ----------------------------- 




## Proc.tiff Function ----------------------------------------------------------

proc.tiff <- function(fit.model,
                      mod.vars,
                      raster.files = NULL,
                      output.options = NULL,
                      factor.levels=NA,
                      make.binary.tif=F,
                      make.p.tif=T,
                      thresh = 0.5,
                      # outfile.p="brt.prob.map.tif",
                      # outfile.bin="brt.bin.map.tif",
                      tsize=2.0,
                      NAval=-3000,
                      fnames=NULL,
                      out,
                      Model){
  
  # vnames,fpath,myfun,make.binary.tif=F,outfile=NA,outfile.bin=NA,output.dir=NA,tsize=10.0,NAval=NA,fnames=NA
  # Written by Alan Swanson, YERC, 6-11-08
  # Revised and Edited by Marian Talbert 2010-2011
  # Description:
  # This function is used to make predictions using a number of .tiff image inputs
  # in cases where memory limitations don't allow the full images to be read in.
  #
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
  # fnames: if the filenames of input files are different from the variable names used in the
  #   prediction model.
  #
  # Modification history:
  # Fixed problem with NA values causing crash 10/2010
  # Included code to produce MESS map and Mod map  8/2011
  # Removed Tiff Directory option as well as some other unused code 10/2011
  #
  # Description:
  # This function reads in a limited number of lines of each image (specified in terms of the
  # size of the temporary predictor dataframe), applies a user-specified
  # prediction function, and stores the results as matrix.  Alternatively, if an
  # output file is specified, a file is written directly to that file in .tif format to
  # the same directory as the input files.  Geographic information from the input images
  # is retained.
  #
  
  # Start of function #
  
  if(is.null(factor.levels)){ factor.levels <- NA }
  # if(is.null(thresh)){ thresh <- 0.5 }
  
  makeMESS <- output.options$MakeMessMap 
  makeMOD <- output.options$MakeModMap
  
  nVars <- length(mod.vars)
  
  if(nVars < 1){
    makeMESS <- FALSE
    message("NOTE: MESS map not generated. Model must have 1 or more predictor variables to make a MESS map.")
  } 
  if(nVars == 1){
    makeMOD <- FALSE 
    message("NOTE: MOD map not generated. Model must have at least 1 predictor variables to make a MOD map.")
  }
  
  ######################################
  # get spatial reference info from existing image file
  
  gi <- GDALinfo(fullnames[1])
  
  dims <- as.vector(gi)[1:2]
  ps <- as.vector(gi)[6:7]
  ll <- as.vector(gi)[4:5]
  pref<-attr(gi,"projection")
  
  RasterInfo=raster(fullnames[1])
  RasterInfo@file@datanotation<-"FLT4S"
  NAval<- -3.399999999999999961272e+38
  
  #To remove use of the Raster package I need to see if rgdal handles area or point correctly
  if(!is.na(match("AREA_OR_POINT=Point",attr(gi,"mdata")))){                                                                          
    xx<-RasterInfo  #this shifts by a half pixel
    nrow(xx) <- nrow(xx) - 1
    ncol(xx) <- ncol(xx) - 1
    rs <- res(xx)
    xmin(RasterInfo) <- xmin(RasterInfo) - 0.5 * rs[1]
    xmax(RasterInfo) <- xmax(RasterInfo) - 0.5 * rs[1]
    ymin(RasterInfo) <- ymin(RasterInfo) + 0.5 * rs[2]
    ymax(RasterInfo) <- ymax(RasterInfo) + 0.5 * rs[2]
  }
  
  # setting tile size
  MB.per.row<-dims[2]*nvars*32/8/1000/1024
  if(MESS) MB.per.row<-MB.per.row*8 #use more blocks for mess
  nrows<-min(round(tsize/MB.per.row),dims[1])
  bs<-c(nrows,dims[2])
  chunksize<-bs[1]*bs[2]
  tr<-blockSize(RasterInfo,chunksize=chunksize)
  
  FactorInd<-which(!is.na(match(mod.vars,names(factor.levels))),arr.ind=TRUE)
  if((nvars-length(FactorInd))==0) MESS<-MOD<-FALSE #turn this off if only one factor column was selected
  
  #for debugging I'm always using multiple cores
  # multCore<-multCore
  if(tr$n<10 | getRversion()<2.14) multCore<-FALSE #turn off multicore in certian circumstances
  
  if(multCore){
    library(parallel)
    #create some temporary folders    
    if(make.p.tif)
      dir.create(file.path(out$input$output.dir,"ProbTiff"))
    outfile.p=file.path(out$input$output.dir,"ProbTiff","_prob_map.tif")
    if(make.binary.tif)                                                                                         
      outfile.bin=dir.create(file.path(out$input$output.dir,"BinTiff"))
    if(MESS) dir.create(file.path(out$input$output.dir,"MESSTiff"))
    if(MOD)  dir.create(file.path(out$input$output.dir,"ModTiff"))
    if(out$input$ResidMaps)
      dir.create(file.path(out$input$output.dir,"ResidTiff"))
    tile.start<-seq(from=1,to=tr$n,by=ceiling(tr$n/(detectCores()-1))) 
    cl <- makeCluster(detectCores()) 
    parLapply(cl,X=tile.start,fun=parRaster,dims=dims,
              tr=tr,MESS=MESS,MOD=MOD,nvars=nvars,fullnames=fullnames,nvars.final=nvars.final,mod.vars=mod.vars,NAval=NAval,factor.levels=factor.levels,
              model=model,Model=Model,pred.fct=pred.fct,make.binary.tif=make.binary.tif,make.p.tif=make.p.tif,RasterInfo=RasterInfo,outfile.p=outfile.p,
              outfile.bin=outfile.bin,thresh=thresh,nToDo= ceiling(tr$n/(detectCores()-1)),ScriptPath=out$input$ScriptPath,
              mod.vars.final.mod=mod.vars.final.mod,train.dat=out$dat$ma$train$dat,residSmooth=out$mods$auc.output$residual.smooth.fct,
              template=out$dat$input$ParcTemplate,maDir=out$input$ma.name)
    stopCluster(cl)
  }  else{  #multicore is slower for small tiffs so we won't do it and the library is not available prior to 2.14
    #also due to multicore multiinstance R issues we're currently only running it on condor or when running synchronously
    parRaster(start.tile=1,dims=dims,
              tr=tr,MESS=MESS,MOD=MOD,nvars=nvars,fullnames=fullnames,nvars.final=nvars.final,mod.vars=mod.vars,NAval=NAval,factor.levels=factor.levels,
              model=model,Model=Model,pred.fct=pred.fct,make.binary.tif=make.binary.tif,make.p.tif=make.p.tif,RasterInfo=RasterInfo,outfile.p=outfile.p,outfile.bin=outfile.bin,thresh=thresh,nToDo=tr$n,ScriptPath=out$       
                input$ScriptPath,mod.vars.final.mod=mod.vars.final.mod,train.dat=out$dat$ma$train$dat,residSmooth=out$mods$auc.output$residual.smooth.fct,
              template=out$dat$input$ParcTemplate,maDir=out$input$ma.name)
  }
  if(length(l<-list.files(dirname(outfile.p),pattern="_prob_map.txt",full.names=TRUE))!=0) unlink(l)
  return(0)
}