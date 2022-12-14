## ------------------------------
## wisdm - apply model functions  
## ApexRMS, April 2022       
## ------------------------------ 


# Calculate MESS Map function ------------------------------------------------
  
  CalcMESS <- function(rast,train.dat){  
    
    # if anything is out of range, return it before calculating sums                       
    min.train <- data.frame(train.dat[1,]) # because we sorted
    max.train <- data.frame(train.dat[nrow(train.dat),])
    output <- data.frame(matrix(data=NA,nrow=nrow(rast),ncol=ncol(rast)))
    names(output) <- names(train.dat)
    for(k in 1:length(min.train)){
      output[,k] <- my.min(rast.val=rast[,k], as.numeric(min.train[,k]), as.numeric(max.train[,k]))
    }
    Indx <- apply(output,1,FUN = function(x){names(which.min(x))})
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
  
## My Min Function -----------------------------------------------------------
  
  my.min <- function(rast.val,min.train,max.train){
    pmin((rast.val-min.train),(max.train-rast.val))/(max.train-min.train)*100
  }

## Vector Sum function -------------------------------------------------------

  vecSum <- function(v,vect){ sum(v > vect) }
    
# Predict Surface function ---------------------------------------------------
  # object <- raster::raster(file.path(ssimTempDir,paste0(modType,"_prob_map.tif")))
  # 
  # Pred.Surface <- function(object, 
  #                          model, 
  #                          filename="", 
  #                          na.rm=TRUE,
  #                          NAval) {
  #   
  #   predrast <- rast(object)
  #   # firstrow <- 1
  #   # firstcol <- 1
  #   # ncols <- ncol(predrast)
  #   # lyrnames <- names(object)
  #   # xylyrnames <- c('x', 'y', lyrnames)
  #   # v <- matrix(NA, ncol=nrow(predrast), nrow=ncol(predrast))
  #   # na.rm <- FALSE
  #   
  #   # tr <- raster::blockSize(predrast, n=raster::nlayers(object)+5)
  #   ablock <- 1:(ncol(object) * nrow(object))
  #   
  #   # napred <- rep(NA, ncol(predrast)*tr$nrows[1])
  #   # predrast <- writeStart(predrast, filename=filename,overwrite=TRUE)
  #   # ############################################################
  #   # for (i in 1:tr$n) {
  #   #   if (i==tr$n) { 
  #   #     ablock <- 1:(ncol(object) * tr$nrows[i])
  #   #     napred <- rep(NA, ncol(predrast) * tr$nrows[i])
  #   #   }
  #   #   rr <- firstrow + tr$row[i] - 1
  #     p <- xyFromCell(predrast, ablock) #, ablock + (tr$row[i]-1) * ncol(predrast)) 
  #     p <- na.omit(p)
  #     blockvals <- data.frame(x=p[,1], y=p[,2])
  #     # if (na.rm) {
  #     #   blockvals <- na.omit(blockvals)		
  #     # }
  #     if (nrow(blockvals) == 0 ) {
  #       predv <- napred
  #     } else {
  #       
  #       #blockvals <- values(object)
  #       predv <- predict(model, blockvals)
  #       predv[is.na(predv)]<-NA
  #     }
  #     # if (na.rm) {  
  #     #   naind <- as.vector(attr(blockvals, "na.action"))
  #     #   if (!is.null(naind)) {
  #     #     p <- napred
  #     #     p[-naind] <- predv
  #     #     predv <- p
  #     #     rm(p)
  #     #   }
  #     # }
  #     
  #     # to change factor to numeric; should keep track of this to return a factor type RasterLayer
  #     predv = as.numeric(predv)
  #     predrast <- writeRaster(predrast, predv, )
  #     #NAvalue(predrast)<-NAval
  #     # print(i)
  #   }
  #   
  #   # predrast <- writeStop(predrast)
  #   
  #   return(predrast)
  # }
  