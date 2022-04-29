## -----------------------------------
## wisdm - data preparation functions  
## ApexRMS, March 2022       
## -----------------------------------

# DATA PREP FUNCTIONS ----------------------------------------------------------

## Split Data for Cross Validation ----------------------------------------------

CrossValidationSplit <- function(data = siteDataWide, # site data
                                 covData = covariatesSheet,
                                 n.folds = 10, # ValidationDataSheet$NumberOfFolds
                                 stratify = FALSE, # ValidationDataSheet$StratifyFolds,
                                 spatSplit = FALSE,
                                 seed = NULL){
  
  #Description:
  #this code takes as input an mds file with the first line being the predictor or
  #response name, the second line an indicator of predictors to include and the
  #third line being paths where tif files can be found.   An output file and a
  #response column must also be specified.  Given a number of folds, a new
  #column is created indicating which fold each observation will be assigned to
  #if a Split Column is also found (test/train split) then only the train portion
  #will be assigned a fold.  Optional stratification by response is also available
  #Background points are ignored
  #by this module (they are read in, written out, but not assigned to cv folds.
  #Output is written to a csv that can be used by the
  #SAHM R modules.
  
  
  if(is.null(seed)) seed <- round(runif(1,min=-((2^32)/2-1),max=((2^32)/2-1)))
  set.seed(seed)
  options(warn=1)
  if(n.folds<=1 | n.folds%%1!=0) stop("n.folds must be an integer greater than 1")
  
  # Read input data and remove any columns to be excluded
  dat <- data
  
  # if there was a test training split this should be used for Evaluation
  # of the final model since cross validation can only be
  # used for model selection
  
  response <- dat$Response
  if(any(response==-9998)) {
    response[response==-9998]<-0
  }
  
  if(sum(as.numeric(response)==0)==0 && !is.null(stratify))
    stop("Cross Validation requires absence data")
  
  # Ignoring background data that might be present in the mds
  
  bg.dat <- dat[response ==-9999,]
  if(dim(bg.dat)[1]!=0){
    dat <- dat[-c(which(response==-9999,arr.ind=TRUE)),]
    response <- response[-c(which(response==-9999,arr.ind=TRUE))]
    bg.dat$Split=""
  }
  
  
  # tagging factors and looking at their levels warning users if their factors have few levels
  factor.cols <- which(names(dat) %in% covData$CovariateName[which(covData$IsCategorical == T)])
  if(length(factor.cols)!=0){
    for (i in 1:length(factor.cols)){
      factor.table <- table(dat[,factor.cols[i]])
      if(any(factor.table<10))
      {warning(paste("Some levels for the categorical predictor ",
                     names(dat)[factor.cols[i]],
                     " do not have at least 10 observations.\n",
                     "you might want to consider removing or reclassifying this predictor before continuing.\n",
                     "Factors with few observations can cause failure in model fitting when the data is split and cannot be reilably used in training a model.",sep=""))
        factor.table <- as.data.frame(factor.table)
        colnames(factor.table)<-c("Factor Name","Factor Count")
        cat(paste("\n",names(dat)[factor.cols[i]],"\n"))
        print(factor.table)
        cat("\n\n")
      }
    }
  }
  # this splits the training set
  if(any(dat$UseInModelEvaluation == TRUE)){ 
    split.mask <- dat$UseInModelEvaluation == FALSE # FALSE = train, TRUE = test 
    index <- seq(1:nrow(dat))[split.mask]
  } else { split.mask <- index <- seq(1:nrow(dat)) }
  
  
  dat[,ncol(dat)+1] <- NA
  
  if(stratify & !spatSplit){ # if stratify = T and spatial split = F
    
    for(i in 1:length(names(table(response)))){
      index.i <- index[response[split.mask] == names(table(response))[i]]
      index.i <- index.i[order(runif(length(index.i)))]
      dat[index.i,ncol(dat)] <- c(rep(seq(1:n.folds), each = floor(length(index.i)/n.folds)),
                                  sample(seq(1:n.folds), size = length(index.i) %% n.folds, replace=FALSE))
    }
    # } else if(spatSplit){ # if spatial split = T
    #   X<-as.numeric(dat$X[split.mask])
    #   Y<-as.numeric(dat$Y[split.mask])
    #   Xbreak<-c(quantile(X,probs=c(.333,.666)),max(X))
    #   Ybreak<-c(quantile(Y,probs=c(.333,.666)),max(Y))
    #   CVsplit<-apply(outer(X,Xbreak,"<"),1,sum)
    #   if(any(CVsplit==0)) CVsplit[CVsplit==0]<-1
    #   CVsplitY<-apply(outer(Y,Ybreak,"<"),1,sum)
    #   if(any(CVsplitY==0)) CVsplitY[CVsplitY==0]<-1
    #   dat[index,ncol(dat)]<-as.numeric(as.factor(paste(CVsplit,CVsplitY)))
  } else { # if stratify = F and spatial split = F
    index <- index[order(runif(length(index)))]
    dat[index,ncol(dat)] <- c(rep(seq(1:n.folds),each = floor(length(index)/n.folds)),
                              sample(seq(1:n.folds), size = length(index) %% n.folds, replace=FALSE))
  }
  names(dat)[ncol(dat)]<-"Split"
  
  if(any(table(dat$Response,dat$Split)<2))
    stop("Some combinations of the response and data split is nearly empty, please use a stratified nonspatial split")
  if(dim(bg.dat)[1] != 0) {
    # names(bg.dat) <- names(dat)
    dat <- rbind(dat, bg.dat)
  }
  
  siteSplits <- dat[,c("SiteID", "Split")]
  
  return(dat)   
  
}
