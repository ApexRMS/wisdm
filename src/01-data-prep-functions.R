## -----------------------------------
## wisdm - data preparation functions  
## ApexRMS, March 2022       
## -----------------------------------

# DATA PREP FUNCTIONS ----------------------------------------------------------

## Split data for Training/Testing ---------------------------------------------

testTrainSplit <- function(inputData,             # dataframe with field data and Covariate Values 
                           trainProp,             # proportion of data to use in model training
                           factorVars = NULL,     # names of categorical variables (if any)
                           ratioPresAbs = NULL){  # ratio of presence to absence points that should be used
                           
  
  # Description:
  # Given a training proportion, each observation is assigned to the
  # test or training split.The split will be balanced with respect to the response 
  # in that the training proportion will be matched as close as possible for each 
  # level of the response.
  # An optional parameter, ratioPresAbs, can be used to specify if there is a certain
  # ratio of presence to absence points that should be used this ensures the given
  # ratio for all data used in both the test and train split. For categorical data,
  # all responses greater than or equal to 1 will be used to meet this ratio.
  # This option reduces the sample size as some data must be thrown away to meet
  # the constraint of having the desired proportion.
  # Background points are ignored by this module (they are read in and written out, 
  # but not assigned to either the test or training split). 
  
  if(trainProp<=0 | trainProp>1){ stop("Training Proportion must be a number between 0 and 1 excluding 0") }
  if(!is.null(ratioPresAbs)) {
    if(ratioPresAbs<=0) stop("The ratio of presence to absence (ratioPresAbs) must be a \ngreater than 0")}
  
  #Read input data and remove any columns to be excluded
  dat <- inputData
  response <- dat$Response
  
  if(sum(as.numeric(response)==0)==0 && !is.null(ratioPresAbs)) stop("The ratio of presence to absence cannot be set with only presence data")
  
  # Ignore background data if present
  bg.dat <- dat[response == -9999,]
  
  if(dim(bg.dat)[1] != 0){
    dat <- dat[-c(which(response == -9999, arr.ind=TRUE)),]
    response <- response[-c(which(response==-9999, arr.ind=TRUE))]
  }
  
  ## Warnings section ##
  
  # tagging factors and looking at their levels warning users if their factors have few levels
  if(!is.null(factorVars)){
    for (i in 1:length(factorVars)){
      factor.table <- table(dat[,factorVars[i]])
      if(any(factor.table<10)){ warning(paste("Some levels for the categorical predictor ",factorVars[i]," do not have at least 10 observations.\n",
                                              "You might want to consider removing or reclassifying this predictor before continuing.\n",
                                              "Factors with few observations can cause failure in model fitting when the data is split and cannot be reilably used in training a model.",sep=""))
        factor.table <- as.data.frame(factor.table)
        colnames(factor.table) <- c("Factor Name","Factor Count")
        cat(paste("\n",factorVars[i],"\n"))
        print(factor.table)
        cat("\n\n")
      }
    }
  }
  if(length(response)<100){ stop("A test training split is not advisable for less than 100 observations.  Consider-cross validation as an alternative.")}
  if(length(response)<200){
    warning(paste("There are less than 200 observations. Cross-validation might be preferable to a test:",
                  "training split \n weigh the decision while keeping in mind the number of predictors being considered: ", ncol(dat)-7,sep=""))
  }
  if(all(na.omit(response) %in% 0:1) & any(table(response)<10)){
    stop("Use of a test training split is not recommended when the dataset contains less than 10 presence or absence points")
  }
  ## End warnings section ##
  
  temp <- if(!is.null(ratioPresAbs)){(sum(response>=1)/sum(response==0) == ratioPresAbs)}
  if(is.null(temp)){ temp <- FALSE }
  
  if(is.null(ratioPresAbs) | temp){
    # Randomly sample presence absence or counts as close to the size of the training proportion as possible
    
    # iterate through each unique response and randomly assign the trainProp to be in the training split
    TrainSplit <- numeric()
    for(i in sort(as.numeric(unique(response)))){
      TrainSplit <- c(TrainSplit, sample(which(response == i, arr.ind=TRUE), size = round(sum(response==i)*trainProp)))
    }
    
    # take everything not in the training set for the test set
    TestSplit <- which(!(seq(1:length(response))) %in% TrainSplit, arr.ind=TRUE)
    
    dat$UseInModelEvaluation[TrainSplit] <- F
    dat$UseInModelEvaluation[TestSplit] <- T

  } else {  # now considering if there is a desired ratio of presence to absence points
    
    if(sum(response>=1)/sum(response==0) >= ratioPresAbs){
      
      # first determine how many presence points to remove
      TotToRmv <- (sum(response>=1)-ratioPresAbs*sum(response==0))
      if(TotToRmv/sum(response>=1) > 0.5){
      warning("******************************************\n**** Over 50% of the 
      presence points were removed to meet the desired ratio of presence to absence
              \n******************************************")}
      
      # determine the number of each count to remove weighted by the response and then remove these
      EachToRmv <- round(TotToRmv * table(response[response!= 0])/sum(response!= 0))
      ByCount <- split(cbind(which(response != 0, arr.ind=TRUE)),f = response[response!=0])
      
      # sampling function ----
      sam <- function(x,size){
        if(length(x)==1 & size==1){ return(x)
        } else { sample(x=x,size=size) }
      }
      
      RmvIndx <- as.vector(unlist(mapply(sam,x=ByCount,size=EachToRmv)))
      KeepIndx <- seq(1:length(response))[-c(RmvIndx)]
      respKeep <- response[KeepIndx]
      names(respKeep) <- KeepIndx
      
      # now break these into a train an test split while
      TrainSplit <- numeric()
      
      for(i in sort(as.numeric(unique(respKeep)))){
        TrainSplit <- c(TrainSplit, sample(names(respKeep[respKeep==i]), size=round(sum(respKeep==i)*trainProp)))
      }
      TrainSplit <- as.numeric(TrainSplit)
      # Take everything not in the training set or in the remove list for the test set
      TestSplit <- seq(from=1,to=length(response))[-c(c(TrainSplit,RmvIndx))]
      
      dat$UseInModelEvaluation[TrainSplit] <- F
      dat$UseInModelEvaluation[TestSplit] <- T
      
    }
    
    if(sum(response>=1)/sum(response==0) < ratioPresAbs){
      
      # first balance all responses greater than 1
      TrainSplit <- numeric()
      for(i in sort(as.numeric(unique(response[response!=0])))){
        TrainSplit <- c(TrainSplit, sample(which(response==i, arr.ind=TRUE), size = round(sum(response==i)*trainProp)))
      }
      if((sum(response==0)-sum(response>=1)/ratioPresAbs)/sum(response==0) > 0.5){
      warning("******************************************\n**** Over 50% of the
              absence points were removed to meet the desired ratio of presence to 
              absence \n******************************************")}
      
      # sample the right number of absence points for the train split
      TrainSplit <- c(TrainSplit, sample(which(response==0, arr.ind=TRUE), size=round(sum(response>=1)*(1/ratioPresAbs)*trainProp)))
      
      # take everything not in the training set for the test set
      TestSplit <- which(!(seq(1:length(response))) %in% TrainSplit, arr.ind=TRUE)
      
      # now sample some points to remove so we have the correct proportion
      temp <- sample(which(TestSplit %in% which(response==0, arr.ind=TRUE), arr.ind=TRUE),
                   size = round(sum(TestSplit %in% which(response==0, arr.ind=TRUE)) - 
                                  (1-trainProp)*sum(response>=1)*(1/ratioPresAbs)))
      TestSplit <- TestSplit[-c(temp)]
      
      dat$UseInModelEvaluation[TrainSplit] <- F
      dat$UseInModelEvaluation[TestSplit] <- T
      
    }
    
    # dat <- dat[c(TrainSplit,TestSplit),]
    
  }
  
  if(dim(bg.dat)[1] != 0) {
    names(bg.dat) <- names(dat)
    dat <- rbind(dat, bg.dat)
    }
  
  return(dat)

}

## Split data for Cross Validation ---------------------------------------------

crossValidationSplit <- function(inputData,         # dataframe with field data and Covariate Values
                                 factorVars,        # names of categorical variables (if any)
                                 nFolds = 10,       # ValidationDataSheet$NumberOfFolds
                                 stratify = FALSE,  # ValidationDataSheet$StratifyFolds,
                                 spatSplit = FALSE){
  
  # Description:
  # this code takes as input an mds file with the first line being the predictor or
  # response name, the second line an indicator of predictors to include and the
  # third line being paths where tif files can be found.   An output file and a
  # response column must also be specified.  Given a number of folds, a new
  # column is created indicating which fold each observation will be assigned to
  # if a Split Column is also found (test/train split) then only the train portion
  # will be assigned a fold.  Optional stratification by response is also available
  # Background points are ignored
  # by this module (they are read in, written out, but not assigned to cv folds.
  
  # cross validation can only be used for model selection
  
  options(warn=1)
  if(nFolds<=1 | nFolds%%1!=0) stop("Number of folds must be an integer greater than 1")
  
  # read input data and remove any columns to be excluded
  dat <- inputData
  response <- inputData$Response
  
  if(sum(as.numeric(response)==0)==0 && !is.null(stratify)){
    stop("Cross Validation requires absence data.")
  }
    
  # ignore background data that might (if present)
  bg.dat <- dat[response ==-9999,]
  
  if(dim(bg.dat)[1]!=0){
    dat <- dat[-c(which(response==-9999, arr.ind=TRUE)),]
    response <- response[-c(which(response==-9999, arr.ind=TRUE))]
  }
  
  
  # tagging factors and looking at their levels warning users if their factors have few levels
  if(!is.null(factorVars)){
    for (i in 1:length(factorVars)){
      factor.table <- table(dat[,factorVars[i]])
      if(any(factor.table<10))
      {warning(paste("Some levels for the categorical predictor ", factorVars[i],
                     " do not have at least 10 observations.\n",
                     "you might want to consider removing or reclassifying this predictor before continuing.\n",
                     "Factors with few observations can cause failure in model fitting when the data is split and cannot be reilably used in training a model.",sep=""))
        factor.table <- as.data.frame(factor.table)
        colnames(factor.table)<-c("Factor Name","Factor Count")
        cat(paste("\n",factorVars[i],"\n"))
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
  
  
  if(stratify & !spatSplit){ # if stratify = T and spatial split = F
    
    for(i in 1:length(names(table(response)))){
      index.i <- index[response[split.mask] == names(table(response))[i]]
      index.i <- index.i[order(runif(length(index.i)))]
      dat$ModelSelectionSplit[index.i] <- c(rep(seq(1:nFolds), each = floor(length(index.i)/nFolds)),
                                            sample(seq(1:nFolds), size = length(index.i) %% nFolds, replace=FALSE))
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
    #   dat$ModelSelectionSplit[index]<-as.numeric(as.factor(paste(CVsplit,CVsplitY)))
  
    } else { # if stratify = F and spatial split = F
    
      index <- index[order(runif(length(index)))]
      dat$ModelSelectionSplit[index] <- c(rep(seq(1:nFolds),each = floor(length(index)/nFolds)),
                                          sample(seq(1:nFolds), size = length(index) %% nFolds, replace=FALSE))
  }
  
  if(all(na.omit(response) %in% 0:1) & any(table(dat$Response,dat$ModelSelectionSplit)<2)){
    stop("Some combinations of the response and data split are nearly empty, please use a stratified nonspatial split")
  }
    
  if(dim(bg.dat)[1] != 0){
    names(bg.dat) <- names(dat)
    dat <- rbind(dat, bg.dat)
  }
  
  return(dat)   
  
}
