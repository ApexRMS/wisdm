## ------------------------------
## wisdm - apply model functions  
## ApexRMS, April 2022       
## ------------------------------ 

# Update breakpoint function ---------------------------------------------------
# Function to time code by returning a clean string of time since this function was last called

updateBreakpoint <- function() { 

  # Calculate time since last breakpoint
  newBreakPoint <- proc.time()
  elapsed <- (newBreakPoint - currentBreakPoint)['elapsed']
  
  # Update current breakpoint
  currentBreakPoint <<- newBreakPoint
  
  # Return cleaned elapsed time
  if (elapsed < 60) {
    return(paste0(round(elapsed), "sec"))
  } else if (elapsed < 60^2) {
    return(paste0(round(elapsed / 60, 1), "min"))
  } else
    return(paste0(round(elapsed / 60 / 60, 1), "hr"))
}


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

  