## ------------------------- 
## sdsim - helper functions  
## ApexRMS, March 2022       
## ------------------------- 

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

# MODEL SELECTION AND VALIDATION FUNCTIONS -------------------------------------

## Predict Function ------------------------------------------------------------

pred.fct <- function(mod,      # mod = the model fit object
                     x,        # x = data to predict for
                     modType){ # modType = one of mars, glm, rf, brt, maxlike at present 
  
  y <- rep(NA,nrow(x))
  
  if(modType %in% c("glm","mars")){
    if("list" %in% class(mod)){
      y <- try(as.vector(predict(mod, x, type="response")),silent=TRUE)
    } else { y <- try(as.vector(predict(mod, x, type="response")),silent=TRUE) }
  } 
  
  # if(modType=="brt"){
  #   # retrieve key items from the global environment #
  #   # make predictions from complete data only #
  #   #y <- rep(NA,nrow(x))
  #   #y[complete.cases(x)] <- predict.gbm(model, x[complete.cases(x),],model$target.trees,type="response")
  #   if(class(model[[1]])=="gbm"){
  #     prd<-function(model,x){
  #       preds <- rep(NA,nrow(x))
  #       preds[complete.cases(x)]<-predict.gbm(model,newdata=x[complete.cases(x),],n.trees=model$target.trees,type="response")
  #       return(preds) 
  #     }         
  #     #getting the predictions from each split of the data then taking out one column and getting the average
  #     lst.preds<-try(lapply(model,FUN=prd,x=x))
  #     y<-try(apply(do.call("rbind",lst.preds),2,mean))
  #   }  else{
  #     # make predictions from full data #
  #     y[complete.cases(x)] <- try(predict.gbm(model,x[complete.cases(x),],model$target.trees,type="response"),silent=TRUE)        
  #   }
  # }
  # if(modType=="rf"){
  #   # retrieve key items from the global environment #
  #   # make predictions from complete data only #
  #   if(class(model[[1]])=="randomForest"){
  #     #getting the predictions from each split of the data then taking out one column and getting the average
  #     lst.preds<-try(lapply(lapply(model,FUN=predict,newdata=x[complete.cases(x),],type="vote"),"[",2))
  #     y[complete.cases(x)]<-try(apply(do.call("rbind",lst.preds),2,mean))
  #     y[y==1]<-max(y[y<1],na.rm=TRUE)
  #     y[y==0]<-min(y[y>0],na.rm=TRUE)
  #   }  else{
  #     y[complete.cases(x)] <- try(as.vector(predict(model,newdata=x[complete.cases(x),],type="vote")[,2]),silent=TRUE)
  #   }  
  # }
  # if(modType=="maxent"){
  #   y[complete.cases(x)]<-try(maxent.predict(model,x[complete.cases(x),]),silent=TRUE)
  # }
  # if(modType=="udc"){
  #   
  #   y[complete.cases(x)]<-try(udc.predict(model,x[complete.cases(x),]),silent=TRUE)
  # }
  # if(class(y)=="try-error") stop("Predicting the response for the new values failed.  One probable cause is that you are trying to predict to factor levels that were not present during model fitting.")
  return(y)
} # end pred.vals function 


## Run Cross Validation --------------------------------------------------------

cv.fct <- function(out,         # out list 
                   nfolds,      # number of cross validation folds
                   sp.no = 1, 
                   prev.stratify = F)
  {
    # function to perform k-fold cross validation
    # with full model perturbation for each subset
    #
    # requires mda library from Cran
    # requires functions mw and calibration
    #
    # takes a mars/glm object produced by mars.glm
    # and first assesses the full model, and then
    # subsets the dataset into nk folds and drops
    # each subset in turn, fitting on remaining data
    # and predicting for withheld data
    #
    # caters for both single species and community models via the argument sp.no
    # for the first, sp.no can be left on its default of 1
    # for community models, sp.no can be varied from 1 to n.spp

    data <- out$data$train  # training data
    n.cases <- nrow(data)
    xdat <- subset(data, select = out$inputVars)
    obs <- data$Response
    preds <- data$predicted
    family <- out$modelFamily
    site.weights <- data$Weight
    
    
    if (family == "binomial" | family == "bernoulli") {
      full.resid.deviance <- calc.deviance(obs, preds, weights = site.weights, family="binomial")
      full.test <- roc(obs, preds)
      full.calib <- calibration(obs, preds)
    }
    
    if (family == "poisson") {
      full.resid.deviance <- calc.deviance(obs, preds, weights = site.weights, family="poisson")                      
      full.test <- cor(obs, preds)
      full.calib <- calibration(obs, preds, family = "poisson")
    }
    
    # set up for results storage
    nk <- nfolds # -1
    subset.test <- rep(0,nk)
    subset.calib <- as.data.frame(matrix(0,ncol=5,nrow=nk))
    names(subset.calib) <- c("intercept","slope","test1","test2","test3")
    subset.resid.deviance <- rep(0,nk)
    
    # now setup for withholding random subsets
    pred.values <- rep(0, n.cases)
    fitted.values <- rep(0, n.cases)
    
    selector <- data$ModelSelectionSplit
    # resp.curves <- vector("list",nk)
    # out.cv <- list()
    # names(resp.curves) <- seq(1:nk)
   
    ## Start cross validation Fold loop 
    
    cor.mat <- matrix(nrow = length(out$inputVars), ncol=nk)
    rownames(cor.mat) <- out$inputVars
    thresh <- NULL
    
    for (i in 1:nk) {
      cat(i," ")
      model.mask <- selector != i  #used to fit model on majority of data
      pred.mask <- selector == i   #used to identify the with-held subset
      # assign("species.subset", obs[model.mask], pos = 1)
      # assign("predictor.subset", xdat[model.mask, ], pos = 1)
      
      # fit new model 
      cv.final.mod <- fitModel(dat = data[model.mask,],
                               out = out,
                               weight = site.weights[model.mask],
                               Fold=i)                       
      
      fitted.values[pred.mask] <- pred.fct(mod = cv.final.mod,
                                           x = xdat[pred.mask,],
                                           modType = out$modType)
      
      
      cor.mat[,i] <- permute.predict(inputVars = out$inputVars, 
                                     dat = xdat[pred.mask,],
                                     obs = obs[pred.mask],
                                     preds = fitted.values[pred.mask],
                                     mod = cv.final.mod,
                                     modType = out$modType)
      
      thresh[i] <- as.numeric(optimal.thresholds(data.frame(ID = 1:nrow(data[model.mask,]),
                                                                         pres.abs = obs[model.mask],
                                                                         pred = pred.fct(mod = cv.final.mod,
                                                                                         x=xdat[model.mask,],
                                                                                         modType = out$modType)),
                                                              opt.methods = out$modOptions$thresholdOptimization))[2]

      y_i <- obs[pred.mask]
      u_i <- fitted.values[pred.mask]
      weights.subset <- site.weights[pred.mask]
      
      if (family == "binomial" | family=="bernoulli") {
        subset.resid.deviance[i] <- calc.deviance(y_i,u_i,weights = weights.subset, family="binomial")
        subset.test[i] <- roc(y_i,u_i)
        subset.calib[i,] <- calibration(y_i, u_i)
      }
      
      if (family=="poisson"){
        subset.resid.deviance[i] <- calc.deviance(y_i,u_i,weights = weights.subset, family="poisson")
        subset.test[i] <- cor(y_i, u_i)
        subset.calib[i,] <- calibration(y_i, u_i, family = family)
      }
      
    } # end of Cross Validation Fold Loop
    
    
    cat("","\n")
    
    u_i <- fitted.values
    
    if (family =="binomial" | family == "bernoulli") {
      cv.resid.deviance <- calc.deviance(obs, u_i, weights = site.weights, family="binomial")
      cv.test <- roc(obs, u_i)
      cv.calib <- calibration(obs, u_i)
    }
    
    if (family == "poisson"){
      cv.resid.deviance <- calc.deviance(obs, u_i, weights = site.weights, family="poisson")
      cv.test <- cor(obs, u_i)
      cv.calib <- calibration(obs, u_i, family = "poisson")
    }
    
    subset.test.mean <- mean(subset.test)
    subset.test.se <- sqrt(var(subset.test))/sqrt(nk)
    
    subset.test <- list(test.scores = subset.test, subset.test.mean = subset.test.mean,
                        subset.test.se = subset.test.se)
    
    subset.calib.mean <- apply(subset.calib[,c(1:2)],2,mean)
    names(subset.calib.mean) <- names(subset.calib)[c(1:2)] # mean only of parameters
    
    subset.calib <- list(subset.calib = subset.calib,
                         subset.calib.mean = subset.calib.mean)
    
    subset.deviance.mean <- mean(subset.resid.deviance)
    subset.deviance.se <- sqrt(var(subset.resid.deviance))/sqrt(nk)
    
    subset.deviance <- list(subset.deviances = subset.resid.deviance, subset.deviance.mean = subset.deviance.mean,
                            subset.deviance.se = subset.deviance.se)
    
    cv.list <- list(full.resid.deviance = full.resid.deviance, full.test = full.test, full.calib = full.calib, 
                  pooled.deviance = cv.resid.deviance, pooled.test = cv.test, pooled.calib = cv.calib,
                  subset.deviance = subset.deviance, subset.test = subset.test, subset.calib = subset.calib,
                  cor.mat = cor.mat, thresh=thresh) # , resp.curves=resp.curves,
    out$cvResults <- cv.list
    return(out)
}


### Calculate Deviance function ------------------------------------------------

calc.deviance <- function(obs,   # observed response
                          preds, # predicted response 
                          weights = rep(1,length(obs)), 
                          family = "binomial", 
                          calc.mean = TRUE)
  {
    # function to calculate deviance given two vectors of raw and fitted values
    # requires a family argument which is set to binomial by default
    
    if (length(obs) != length(preds)){ stop("observations and predictions must be of equal length") }
      
    if (family == "binomial" | family == "bernoulli") {
      # preds[preds==0] <- 0.000000001
      deviance.contribs <- (obs * log(preds)) + ((1-obs) * log(1 - preds))
      deviance <- -2 * sum(deviance.contribs * weights, na.rm = T)
    }
    
    if (family == "poisson" | family == "Poisson") {
      deviance.contribs <- ifelse(obs == 0, 0, (obs * log(obs/preds))) - (obs - preds)
      # deviance.contribs[which(deviance.contribs==Inf)] <- 1
      deviance <- 2 * sum(deviance.contribs * weights)
    }
    
    if (family == "laplace") { deviance <- sum(abs(obs - preds)) }
    
    if (family == "gaussian") { deviance <- sum((obs - preds) * (obs - preds)) }
    
    if (calc.mean) deviance <- deviance/length(obs)
    
    return(deviance)
    
  }

### ROC Function ---------------------------------------------------------------

roc <- function (obs,    # observed response
                 preds)  # predicted response
                 
  {
    # code adapted from Ferrier, Pearce and Watson's code
    # see:
    # Hanley, J.A. & McNeil, B.J. (1982) The meaning and use of the area
    # under a Receiver Operating Characteristic (ROC) curve.
    # Radiology, 143, 29-36
    #
    # Pearce, J. & Ferrier, S. (2000) Evaluating the predictive performance
    # of habitat models developed using logistic regression.
    # Ecological Modelling, 133, 225-245.
    # this is the non-parametric calculation for area under the ROC curve,
    # using the fact that a MannWhitney U statistic is closely related to
    # the area
    
    if (length(obs) != length(preds))
      stop("obs and preds must be equal lengths")
    n.x <- length(obs[obs == 0])
    n.y <- length(obs[obs == 1])
    xy <- c(preds[obs == 0], preds[obs == 1])
    rnk <- rank(xy)
    wilc <- ((n.x * n.y) + ((n.x * (n.x + 1))/2) - sum(rnk[1:n.x]))/(n.x *
                                                                       n.y)
    return(round(wilc, 4))
  }

### Calibration function -------------------------------------------------------

calibration <- function(obs,   # observed response
                        preds, # predicted response
                        family = "binomial")
  {
    # calculates calibration statistics for either binomial or count data
    # but the family argument must be specified for the latter
    # a conditional test for the latter will catch most failures to specify
    # the family

    if (family == "bernoulli"){ family <- "binomial"}
      pred.range <- max(preds) - min(preds)
    
    if(pred.range > 1.2 & family == "binomial") {
      print(paste("range of response variable is ", round(pred.range, 2)), sep = "", quote = F)
      print("check family specification", quote = F)
      return()
    }
    if(family == "binomial") {
      pred <- preds + 1e-005
      pred[pred >= 1] <- 0.99999
      mod <- glm(obs ~ log((pred)/(1 - (pred))), family = binomial)
      lp <- log((pred)/(1 - (pred)))
      a0b1 <- glm(obs ~ offset(lp) - 1, family = binomial)
      miller1 <- 1 - pchisq(a0b1$deviance - mod$deviance, 2)
      ab1 <- glm(obs ~ offset(lp), family = binomial)
      miller2 <- 1 - pchisq(a0b1$deviance - ab1$deviance, 1)
      miller3 <- 1 - pchisq(ab1$deviance - mod$deviance, 1)
    }
    if(family == "poisson") {
      mod <- glm(obs ~ log(preds), family = poisson)
      lp <- log(preds)
      a0b1 <- glm(obs ~ offset(lp) - 1, family = poisson)
      miller1 <- 1 - pchisq(a0b1$deviance - mod$deviance, 2)
      ab1 <- glm(obs ~ offset(lp), family = poisson)
      miller2 <- 1 - pchisq(a0b1$deviance - ab1$deviance, 1)
      miller3 <- 1 - pchisq(ab1$deviance - mod$deviance, 1)
    }
    calibration.result <- c(mod$coef, miller1, miller2, miller3)
    names(calibration.result) <- c("intercept", "slope", "testa0b1", "testa0|b1", "testb1|a")
    return(calibration.result)
}

### Permute Predict function ---------------------------------------------------

permute.predict <- function(inputVars, # input variables for model fitting
                            dat,       # variable site data 
                            obs,       # response (observed) values
                            preds,     # predicted values
                            mod,       # fit model 
                            modType    # model type
                            ){

  AUC <- matrix(NA, nrow = length(inputVars),ncol=5) 
  
  for(j in 1:5){ # do the permutation 5 times to remove some of the random chance component
    for (i in 1:length(inputVars)){
      indx <- match(inputVars[i],names(dat))
      dat.i <- dat
      dat.i[,indx] <- dat.i[sample(1:dim(dat)[1]),indx]
      options(warn=-1)
      new.pred <- as.vector(pred.fct(mod = mod, x = dat.i, modType = modType))
      # have to use ROC here because auc in presence absence incorrectly assumes auc will be greater than .5
      AUC[i,j] <- roc(obs,new.pred)
      options(warn=0)
    } 
  }
  AUC <- apply(AUC,1,mean)
  return(AUC)
}

# MODEL OUTPUT FUNCTIONS -------------------------------------------------------

## Make AUC function -----------------------------------------------------------

make.auc.plot.jpg <- function(out = out){
  
  standResidualFile <- paste0(out$tempDir, out$modType, "_StandardResidualPlots.png")
  variableImportanceFile <- paste0(out$tempDir, out$modType, "_VariableImportance.png")
  confusionMatrixFile <- paste0(out$tempDir, out$modType, "_ConfusionMatrix.png")
  ROCAUCFile <- paste0(out$tempDir, out$modType, "_ROCAUCPlot.png")  
  AUCPRFile <- paste0(out$tempDir, out$modType, "_AUCPRPlot.png")
  calibrationFile <- paste0(out$tempDir, out$modType, "_CalibrationPlot.png")
  residualPlotFile <- paste0(out$tempDir, out$modType, "_PoissonResidulePlot.png")
  
  if(!out$validationOptions$SplitData & !out$validationOptions$CrossValidate){ splitType <- "none" }
  if(out$validationOptions$SplitData & !out$validationOptions$CrossValidate){ splitType <- "train/test" }
  if(out$validationOptions$SplitData & out$validationOptions$CrossValidate){ splitType <- "cv/train/test"}
  if(!out$validationOptions$SplitData & out$validationOptions$CrossValidate){ splitType <- "cv"}

  ### Calculate threshold on train/testing splits ###
  
  if(out$modelFamily != "poisson"){
    
    # To Do: link opt.method to defined Threshold Optimization Method in UI - currently set to default
    out$trainThresh <- as.numeric(optimal.thresholds(data.frame(ID = 1:length(out$data$train$Response),
                                                                 pres.abs = out$data$train$Response,
                                                                 pred = out$data$train$predicted),
                                                        opt.methods = out$modOptions$thresholdOptimization))[2]
    if(!is.null(out$data$test)){
      
      out$testThresh <- as.numeric(optimal.thresholds(data.frame(ID = 1:length(out$data$test$Response),
                                                                     pres.abs = out$data$test$Response,
                                                                     pred = out$data$test$predicted),
                                                          opt.methods = out$modOptions$thresholdOptimization))[2]
    }
  } else {
    
    out$trainThresh <- NULL
    
  }
  
  ### Create standard residual analysis plots for glm  ###
  
  if(out$modType %in% c("glm","mars") & is.null(out$data$test) & !(out$modType == "mars" & out$pseudoAbs)){
    
    png(standResidualFile, height = 1000, width = 1000)
    par(mfrow = c(2, 2))
    if(out$modType == "glm"){ plot(out$finalMod, cex = 1.5, lwd = 1.5, cex.main = 1.5, cex.lab = 1.5) }
    # if(out$input$script.name == "mars") plot(out$finalMod$glm.list[[1]], cex = 1.5, lwd = 1.5, cex.main = 1.5, cex.lab = 1.5)
    par(mfrow = c(1,1))
    graphics.off()
    
  }
  
  ### Calculate all statistics on test\train or train\cv splits  ###
  
  out$hasSplit <- hasSplit <- (out$pseudoAbs & !out$modType %in% c("glm", "maxent"))
  if(out$validationOptions$CrossValidate){
    Stats <- lapply(out$data$cvSplits$test, calcStat, family = out$modelFamily, has.split = hasSplit)
  } else { Stats <- list() }
  Stats$train <- calcStat(x = out$data$train, family = out$modelFamily, has.split = hasSplit)
  if(out$validationOptions$SplitData){
    Stats$test <- calcStat(x = out$data$test, family = out$modelFamily, has.split = hasSplit)
  }
  
  lst <- Stats
  
  ### Variable Importance Plot ###
  
  if(length(out$inputVars) > 1 & out$modelFamily != "poisson"){
    
    png(variableImportanceFile, height = 1000, width = 1000, pointsize = 13)
    VariableImportance(out = out, auc = lapply(Stats, "[",9)) 
    graphics.off()
    
  }    
  
  ### Confusion Matrix Plot ###
  
  
  if(out$modelFamily != "poisson"){
    
    png(file = confusionMatrixFile, width = 1000, height = 1000, pointsize = 13)
    confusion.matrix(Stats, out)
    graphics.off()
    
  }
  
  # ### Create residual surface of input data  # TO DO: move to map output section 
  # 
  # if(out$input$ResidMaps){
  #   
  #   if(out$dat$split.label != "eval"){
  #     
  #     residual.smooth.fct <- resid.image(calc.dev(inlst$train$dat$response, inlst$train$pred, inlst$train$weight, family = out$input$model.family)$dev.cont, inlst$train$pred,
  #                                        inlst$train$dat$response, inlst$train$XY$X, inlst$train$XY$Y, inlst$train$weight, out$input$model.family, out$input$output.dir, label = out$dat$split.label, out)
  #   } else {
  #     
  #     residual.smooth.fct <- resid.image(calc.dev(inlst$test$dat$response, inlst$test$pred, inlst$test$weight, family = out$input$model.family)$dev.cont, inlst$test$pred,
  #                                        inlst$test$dat$response, inlst$test$XY$X, inlst$test$XY$Y, inlst$test$weight, out$input$model.family, out$input$output.dir, label = out$dat$split.label, out)
  #   }
  # } else {
  #   
  #   residual.smooth.fct=NULL
  #   
  # }
  
  ### AUC and Calibration plot for binomial data ### 
  
  
  if(out$modelFamily %in% c("binomial", "bernoulli")){
    
    png(file = ROCAUCFile , height = 1000, width = 1000, pointsize = 20)
    
    #### ROC AUC Plots ####
    
    TestTrainRocPlot(dat = Stats$train$auc.data, 
                     opt.thresholds = out$trainThresh, 
                     add.legend = FALSE, lwd = 2)
    
    if(splitType == "none"){
      
      legend(x = .8, y = .15, paste("AUC=", round(Stats$train$auc.fit, digits = 3), sep = ""))
    
      } else {
      
      if(splitType == "train/test"){ 
        
        lst <- Stats[c("train","test")]
        
        TestTrainRocPlot(do.call("rbind", lapply(lst, function(lst){lst$auc.data})), 
                         add.roc = TRUE, line.type = 2, color = "red", add.legend = FALSE)
        
        legend(x = .46, y = .24, 
               c(paste("Training Split (AUC=", round(Stats$train$auc.fit, digits = 3), ")", sep = ""),
                 paste("Testing Split  (AUC=", round(Stats$test$auc.fit, digits = 3), ")", sep = "")),
               lty = 2, col = c("black", "red"), lwd = 2, cex = 1.3)
      }
      
      if(splitType %in% c("cv", "cv/train/test")){
        
        lst <- Stats[1:out$validationOptions$NumberOfFolds]
        
        ROC.list <- list(predictions = lapply(lst, function(lst){lst$auc.data$pred}), 
                         labels = lapply(lst, function(lst){lst$auc.data$pres.abs}))
        pred <- prediction(ROC.list$predictions, ROC.list$labels)
        perf <- performance(pred, "tpr", "fpr")
        
        plot(perf, col = "grey82", lty = 3, xlab = "1-Specificity (False Positive)", ylab = "Sensitivity (True Positive)",
             main = "ROC Plot for Cross-Validation", cex.main = 2, cex.axis = 1.4, cex.lab = 1.5)
        
        plot(perf, lwd = 1, avg = "vertical", spread.estimate = "boxplot", add = TRUE)
        
        TestTrainRocPlot(dat = Stats$train$auc.data, 
                         opt.thresholds = out$trainThresh,
                         add.legend = FALSE, lwd = 1.5, add.roc = TRUE, 
                         line.type = 1, col = "red", legend.cex = 2)
        
        points(1-Stats$train$Specf, Stats$train$Sens, pch = 21, cex = 2.5, bg = "red")
        segments(x0 = 0, y0 = 0, x1 = 1, y1 = 1, col = "blue")
        text(x = (.96-Stats$train$Specf), y = Stats$train$Sens + .03, label = round(Stats$train$thresh, digits = 2))
        
        legend(x = .5,y = .24,
               c(paste("Training Split (AUC=", round(Stats$train$auc.fit, digits=3), ")", sep = ""),
                 paste("Cross Validation Mean \n (AUC=", round(mean(unlist(lapply(lst, function(lst){lst$auc.fit}))), digits=3), ")", sep="")), 
               lwd = c(4, 1), lty = c(1, 1), col = c("red", "black"), cex = 1.3)
      }
    }
    
    graphics.off()
    
    #### AUCPR Plots #### 
    
    aucpr.lst <- list()
    mean.vec <- c()
    df <- data.frame()
    
    pl <- ggplot(df) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
            legend.position = "bottom", panel.border = element_rect(fill = 'transparent', color = 'black'), 
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"), 
            axis.title = element_text(size = 15),
            axis.title.x = element_text(vjust = -2),
            plot.title = element_text(hjust = 0.5, face = 'bold')) +
      ylab(expression(atop("Precision", atop(scriptscriptstyle(frac("True Positives", "True Positives +  False Positives")))))) +
      xlab(expression(atop("Recall", atop(scriptscriptstyle(frac("True Positives", "True Positives +  False Negatives")))))) +
      ggtitle('PR Plot for Cross Validation') +
      scale_x_continuous(breaks = seq(0, 1, by = 0.2)) +
      scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
      expand_limits(x = 0, y = 0) 
    
    for(n in 1:(length(lst))){
      tmp <- pr.curve(scores.class0 = lst[[n]]$auc.data$pred, 
                      weights.class0 = lst[[n]]$auc.data$pres.abs, curve = T)
      tmp$curve <- as.data.frame(tmp$curve)
      colnames(tmp$curve) <- c('recall', 'precision', 'thresh')
      aucpr.lst[[paste0(n)]] <- tmp
      mean.vec <- c(mean.vec, tmp$auc.davis.goadrich)
      pl <- pl + geom_line(data = tmp$curve, aes(x = recall, y = precision), 
                           color = 'grey', linetype = 'dashed', alpha = 0.3)
    }
    
    # Generate CV mean #
    cv.df <- data.frame()
    for(i in 1:length(aucpr.lst)){
      tmp <- aucpr.lst[[i]]$curve %>%
      summarise_all(.funs = mean)
      cv.df <- bind_rows(cv.df, tmp)
    }
    
    # Add training data # 
    pr <- pr.curve(scores.class0 = Stats$train$auc.data$pred, 
                   weights.class0 = Stats$train$auc.data$pres.abs, curve = TRUE)
    pr$curve <- as.data.frame(pr$curve)
    colnames(pr$curve) <- c('recall', 'precision', 'thresh')
    
    pl <- pl + suppressWarnings(geom_line(data = as.data.frame(spline(pr$curve$recall, pr$curve$precision)),
                                          aes(x = x, y = y, color = 'train'), size = 1)) +
      suppressWarnings(geom_line(data = as.data.frame(spline(pr$curve$recall, pr$curve$precision)), 
                                 aes(x = x, y = y, color = 'test'), size = 1)) +
      scale_color_manual(labels = c(paste0("Training Split (AUC=", round(Stats$train$auc.pr,3),")"),
                                    paste0("Cross Validation Mean (AUC=", mean(round(mean.vec, digits = 3)), ")", sep='')), 
                         values=c('red', 'grey')) +
      theme(legend.title = element_blank(), legend.box.background = element_rect(colour = "black"),
            legend.position = c(0.6, 0.15), legend.key = element_rect(fill = NA))
    
    ggsave(filename = AUCPRFile, 
           plot = pl, width = 12.5, height = 12.5, units = 'cm', dpi = 300)
    
   # TO DO: should above code should only be created if CV ^^^ put into if statement??
    
    #### Calibration Plot #### 
    
    png(file = calibrationFile, height = 1000, width = 1000, pointsize = 20)
    cal.results <- switch(splitType, 
                          "none" = Stats$train$calibration.stats, 
                          "train/test" = Stats$test$calibration.stats, 
                          "cv" = apply(do.call("rbind", lapply(lst, function(lst){lst$calibration.stats})), 2, mean),
                          "cv/train/test" = apply(do.call("rbind", lapply(lst, function(lst){lst$calibration.stats})), 2, mean))
    
    # Calibration plot (this often gives warnings about probabilities numerically 0 or 1)
    a <- do.call("rbind", lapply(lst, function(lst){lst$auc.data}))
    if(out$pseudoAbs){
      if(!(out$modType %in% c("glm"))){
        absn <- which(a$pres.abs == 0, arr.ind = TRUE)
        samp <- sample(absn, size = min(table(a$pres.abs)), replace = FALSE) 
      }
      
      p.plt <- try(pocplot(a$pred[a$pres.abs == 1], a$pred[a$pres.abs == 0],
                           title = paste0("Presence Only Calibration Plot for \n", 
                                         switch(splitType, 
                                                "none" = "Training Data", 
                                                "train/test" = "Test Split",
                                                "cv" = "Cross Validation Split",
                                                "cv/train/test" = "Cross Validation Split"))), silent = TRUE)
      if(class(p.plt) == "try-error"){
        
        par(mfrow = c(2, 1))
        hist(a$pred[a$pres.abs == 1], freq = TRUE, col = "red", 
             xlim = range(a$pred), xlab = "Predicted Probability", main = "Presence")
        hist(a$pred[a$pres.abs == 0], freq = TRUE, col = "blue", 
             xlim = range(a$pred), xlab = "Predicted Probability", main = "Available")
      }
    } else {
      
      pacplot(pred = a$pred, pa = a$pres.abs, 
              title = paste0("Calibration Plot for ",
                             switch(splitType, 
                                    "none" = "Training Data", 
                                    "train/test" = "Test Split",
                                    "cv" = "Cross Validation Split",
                                    "cv/train/test" = "Cross Validation Split")))
    }
  dev.off()
  }
  
  ### Residual plots for poisson data ###
  
  if(out$modelFamily %in% c("poisson")){
    
    png(file = residualPlotFile)
    par(mfrow = c(2, 2))
    
    plot(log(Stats$train$auc.data$pred[Stats$train$auc.data$pred != 0]),
         (Stats$train$auc.data$pres.abs[Stats$train$auc.data$pred!=0]-Stats$train$auc.data$pred[Stats$train$auc.data$pred!=0]),
         xlab="Predicted Values (log scale)", ylab = "Residuals", main = "Residuals vs Fitted", ylim = c(-3, 3))
    abline(h = 0, lty = 2)
    
    panel.smooth(log(Stats$train$auc.data$pred[Stats$train$auc.data$pred != 0]),
                 (Stats$train$auc.data$pres.abs[Stats$train$auc.data$pred!=0]-Stats$train$auc.data$pred[Stats$train$auc.data$pred!=0]))
    
    if(out$modType != "rf"){
      
      # this is the residual plot from glm but I don't think it will work for anything else
      qqnorm(residuals(out$finalMod), ylab = "Std. deviance residuals")
      qqline(residuals(out$finalMod))
      yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name("Std. Deviance Resid"))))
      
      plot(log(Stats$train$auc.data$pred[Stats$train$auc.data$pred != 0]),
           sqrt((abs(residuals(out$finalMod, 
                               type = "deviance")[Stats$train$auc.data$pred != 0]))),
           xlab = "Predicted Values (log Scale)", ylab = yl)
    }
    graphics.off()
  }
  
  ### Text Output ###
  
  # capture.output(cat("\n\n============================================================",
  #                    "\n\nEvaluation Statistics"), file = paste0(out$tempDir, out$modType, "_output.txt"), append = TRUE)
  # 
  # train.stats = list()
  # 
  # if(splitType == "none"){ train.stats <- Stats
  # } else {
  #   train.stats$train=Stats$train
  # }
  # 
  # capture.stats(train.stats, file.name =paste0(out$tempDir, out$modType, "_output.txt"), 
  #               label = "train", family = out$modelFamily, opt.methods = out$modOptions$thresholdOptimization, out)
  # 
  # if(out$dat$split.type != "none"){
  #   
    capture.output(cat("\n\n============================================================",
                       "\n\nEvaluation Statistics"), file = paste0(out$tempDir, out$modType, "_output.txt"), append = TRUE)
    
    capture.stats(Stats, file.name = paste0(out$tempDir, out$modType, "_output.txt"), 
                  label = splitType, family = out$modelFamily, opt.methods = out$modOptions$thresholdOptimization, out)
    
  # }
  
    
    ## STOPPED HERE ------------------------------------------------------------
  ### getting statistics along with appropriate names into a data frame for creating the appended output
  parent <- dirname(out$input$output.dir)
  
  if(out$input$model.family %in% c("binomial", "bernoulli")){
    
    csv.stats <- lapply(Stats, function(lst){return(c("", "", lst$correlation, lst$pct.dev.exp,lst$Pcc, lst$auc.fit, lst$auc.pr,lst$Tss))})
    stat.names <- c("Correlation Coefficient", "Percent Deviance Explained", "Percent Correctly Classified", "AUC", "AUCPR", "True Skill Stat")
  } else {
    
    csv.stats <- lapply(Stats, function(lst){return(c("", "", lst$correlation, lst$pct.dev.exp, lst$prediction.error/100))})
    stat.names <- c("Correlation Coefficient", "Percent Deviance Explained", "Prediction Error")
    
  }
  
  csv.vect <- c(t(t(as.vector(unlist(csv.stats[train.mask])))), if(out$dat$split.type != "none") unlist(csv.stats[-c(train.mask)]))
  csv.vect[seq(from = 2, by = length(csv.vect)/length(Stats), length = length(Stats))] <- if(out$dat$split.type == "none"){
    
    "Train"
  } else {
    
    c("Train",names(lst))
    
  }
  x = data.frame(cbind(rep(c("", "", stat.names), times = length(Stats)), csv.vect), row.names = NULL)
  
  Header<-cbind(c("", "Original Field Data", "Field Data Template", "PARC Output Folder", "PARC Template", "Covariate Selection Name", ""),
                c(basename(out$input$output.dir), out$dat$input$OrigFieldData, out$dat$input$FieldDataTemp, out$dat$input$ParcOutputFolder,
                  basename(out$dat$input$ParcTemplate), ifelse(length(out$dat$input$CovSelectName) == 0,"NONE", out$dat$input$CovSelectName), ""))
  
  AppendOut(compile.out = out$input$Append.Dir, Header, x, out, Parm.Len = length(stat.names), parent = parent, split.type = out$dat$split.type)
  
  return(list(thresh = train.stats$train$thresh, residual.smooth.fct = residual.smooth.fct))
  
}

### Calculate statistics function ----------------------------------------------

calcStat <- function(x,       # x <- out$data[[i]]
                     family,
                     has.split){
  
  auc.data <- data.frame(ID=1:nrow(x),pres.abs=x$Response, pred=x$predicted)
  
  p.bar <- sum(auc.data$pres.abs * x$Weight) / sum(x$Weight)
  n.pres <- sum(auc.data$pres.abs>=1)
  n.abs <- nrow(auc.data)-n.pres
  # if(has.split & !is.null(x$Split)){ #we're in the train split here for pseudoAbs
  #   dev.vect<-vector()
  #   null.dev<-vector()
  #   for(i in 1:(length(unique(x$Split))-1)){
  #     dev.vect[i]<-calc.deviance(c(rep(0,times=sum(x$Split==i)),x$resp[x$resp==1]),
  #                                c(x$pred[x$Split==i],x$pred[x$resp==1]))
  #     null.dev[i]<-calc.deviance(c(rep(0,times=sum(x$Split==i)),x$resp[x$resp==1]),
  #                                rep(mean(c(rep(1,times=sum(x$resp==1)),rep(0,times=sum(x$Split==i)))),times=(sum(x$resp==1)+sum(x$Split==i))))
  #   }
  #   null.dev<-mean(null.dev)
  #   dev.fit<-mean(dev.vect)
  # }
  # else if(has.split & is.null(x$Split)){ # in the test split for pseudoAbs deviance shouldn't be calculated because calibration is wrong
  #   null.dev = NA
  #   dev.fit = NA
  # }
  # else{     
    null.dev <- calc.deviance(auc.data$pres.abs, rep(p.bar,times=length(auc.data$pres.abs)), x$Weight, family=family) #*nrow(x$dat)
    if(is.nan(null.dev)){null.dev <- NA}
    
    dev.fit <- calc.deviance(auc.data$pres.abs, auc.data$pred, x$Weight, family=family) #*nrow(x$dat) Elith does not include this it might cause a weighting issue when averaging I'm not sure if I should include it
    if(is.nan(dev.fit)){dev.fit <- NA}
  # }
  
  dev.exp <- null.dev - dev.fit
  pct.dev.exp <- dev.exp/null.dev*100
  correlation <- cor(auc.data$pres.abs, auc.data$pred)
  
  # have to use roc here because auc in the PresenceAbsence package incorrectly assumes that the value must be greater than .5
  # this isn't necessarily true for an independent evaluation set
  
  auc.fit <- roc(auc.data$pres.abs, auc.data$pred)
  auc.pr <- pr.curve(scores.class0 = auc.data$pred, weights.class0 = auc.data$pres.abs)[[3]]
  
  calibration.stats <- calibration(auc.data$pres.abs, auc.data$pred, family =family)
  
  if(family %in% c("binomial","bernoulli")){
    cmx <- cmx(auc.data,threshold=out$trainThresh)
    PCC <- pcc(cmx,st.dev=F)*100
    SENS <- sensitivity(cmx,st.dev=F)
    SPEC <- specificity(cmx,st.dev=F)
    KAPPA <- Kappa(cmx,st.dev=F)
    TSS <- SENS+SPEC-1
    
    return(list(n.pres=n.pres, n.abs=n.abs,
                null.dev=null.dev, dev.fit=dev.fit, dev.exp=dev.exp,
                pct.dev.exp=pct.dev.exp,
                correlation=correlation,
                auc.data=auc.data, auc.fit=auc.fit, auc.pr=auc.pr,
                Cmx=cmx, Pcc=PCC, Sens=SENS, Specf=SPEC, Kappa=KAPPA, Tss=TSS,
                calibration.stats=calibration.stats, thresh=out$trainThresh))
  }
  
  if(family == "poisson"){
    prediction.error <- sum((auc.data$pres.abs-auc.data$pred)^2)
    return(list(n.pres=n.pres,n.abs=n.abs,
                null.dev=null.dev, dev.fit=dev.fit, dev.exp=dev.exp, 
                pct.dev.exp=pct.dev.exp, 
                correlation=correlation,
                auc.data=auc.data, auc.fit=auc.fit, auc.pr=auc.pr,
                prediction.error=prediction.error, calibration.stats=calibration.stats))
  }
  
}

### Variable Importance function -----------------------------------------------

VariableImportance <- function(out,  # out list
                               auc)  # list of auc.fit values output by calcStat function
  { 
  
  # this relatively complicated structure is used for cross validation variable importance plots
  # This function can produce some strange results for random Forest if it is seriously overparameterized
  # in which case permutation in any single value might not lead to a drop in auc so all variables have
  # approximately zero importance all other models drop variables and haven't shown this
  
  cor.mat <- matrix(nrow = length(out$inputVars), ncol = length(auc))
  
  # # remove the response column
  # if(out$modType == "rf") {
  #   trainPred <- pred.fct(model=out$mods$final,x=out$dat$ma$train$dat,Model=out$input$script.name)
  #   auc$train<-roc(out$dat$ma$train$dat[,1],trainPred)
  # } else 
  
  # add 1 to the drop in AUC to ensure it's greater than zero then I can normalize
  
  cor.mat[,1] <- unlist(auc$train) - permute.predict(inputVars = out$inputVars,
                                                              dat = out$data$train, 
                                                              preds = out$data$train$predicted,
                                                              obs = out$data$train$Response,
                                                              mod = out$finalMod, 
                                                              modType = out$modType)
  cnames <- "train"
  
  # cor.mat[,ncol(cor.mat)]<-cor.mat[,ncol(cor.mat)]/sum(cor.mat[,ncol(cor.mat)])
  
  if(out$validationOptions$SplitData){
    cor.mat[,2] <- unlist(auc$test) - permute.predict(inputVars = out$inputVars,
                                                      dat = out$data$test,
                                                      preds = out$data$test$predicted,
                                                      obs = out$data$test$Response,
                                                      mod = out$finalMod,
                                                      modType = out$modType)
    cnames <- c(cnames,"test")  
  }
  # variable importance values are not normalized - values provided so people can do what they want
  # cor.mat[,1]<-cor.mat[,1]/sum(cor.mat[,1])
  
  if(out$validationOptions$CrossValidate){
    cor.mat[,(length(cnames)+1):(ncol(cor.mat))] <- -t(apply(out$cvResults$cor.mat,1,"-",as.vector(unlist(auc)[1:out$validationOptions$NumberOfFolds])))
    cnames <- c(cnames,1:out$validationOptions$NumberOfFolds)
  }
    
  
  colnames(cor.mat) <- cnames
  rownames(cor.mat) <- out$inputVars
  
  # TO DO: writing output should be formalized more if it is kept
  variable.importance.csv <- file.path(paste0(ssimTempDir,"\\Data\\",out$modType,"_VariableImportance.csv"))  
  write.table(cbind(predictor=out$inputVars,cor.mat),
              file = variable.importance.csv,
              row.names=FALSE,col.names=TRUE,quote=FALSE,sep=",")
  
  
  # if cross validation we need to avg across folds otherwise plot for each
  # order by the best in the train split
  
  xright <- as.matrix(cor.mat[order(cor.mat[,"train"],decreasing=FALSE),])
  ymiddle <- seq(from=0,to=length(out$inputVars),length=nrow(xright))
  offSet <- 0.5 
  
  par(mar=c(6,17,6,0))
  
  if(!out$validationOptions$SplitData & !out$validationOptions$CrossValidate){
    plot(x = c(min(0,min(xright[,"train"])),(max(xright[,"train"])+0.1)),
         y = c(-0.5,(length(out$inputVars)+0.5)),
         type="n",xlab="Importance", ylab="", yaxt="n",
         main="Importance using the change in AUC\nwhen each predictor is permuted",
         cex.lab=3, cex.main=3, cex.axis=2)
    grid()
    rect(xleft=0, ybottom=ymiddle, xright=xright[,"train"], ytop=ymiddle+offSet, col="blue", lwd=2)
    legend("bottomright", legend="Train", fill="blue", bg="white", cex=2.5)
  } 
  
  if(out$validationOptions$SplitData & !out$validationOptions$CrossValidate){
    plot(x = c(min(0,min(xright[,"train"])),(max(xright[,"train"])+0.1)),
         y = c(-0.5,(length(out$inputVars)+0.5)),
         type="n",xlab="Importance", ylab="", yaxt="n",
         main="Importance using the change in AUC\nwhen each predictor is permuted",
         cex.lab=3, cex.main=3, cex.axis=2)
    grid()
    rect(xleft=0, ybottom=ymiddle, xright=xright[,"train"], ytop=ymiddle+offSet, col="blue", lwd=2)
    rect(xleft=0, ybottom=ymiddle-offSet, xright=xright[,"test"], ytop=ymiddle, col="lightblue", lwd=2)
    legend("bottomright", legend=c("Train", "Test"), fill=c("blue","lightblue"), bg="white", cex=2.5)
  } 
  
  if(!out$validationOptions$SplitData & out$validationOptions$CrossValidate){ 
    plot(x = c(min(0,min(xright)),(max(xright)+0.1)),
       y = c(-0.5,(length(out$inputVars)+0.5)),
       type="n",xlab="Importance", ylab="", yaxt="n",
       main="Importance using the change in AUC\nwhen each predictor is permuted",
       cex.lab=3, cex.main=3, cex.axis=2)
    grid()  
    boxplot(t(xright[,-1]), horizontal=TRUE, add=TRUE, at=ymiddle, yaxt="n", col="lightblue", xaxt="n")
    points(y=ymiddle, x=xright[,"train"], cex=3, pch=8, lwd=3, col="darkslateblue")
    legend(x="bottomright", legend=c("CV","Train"),pch=c(22,8),pt.cex=c(3,3.5),pt.lwd=c(2,3),pt.bg=c("lightblue","darkslateblue"),col=c("black","darkslateblue"),cex=2.5)
  }
  
  if(out$validationOptions$SplitData & out$validationOptions$CrossValidate){ 
    plot(x = c(min(0,min(xright)),(max(xright)+0.1)),
         y = c(-0.5,(length(out$inputVars)+0.5)),
         type="n",xlab="Importance", ylab="", yaxt="n",
         main="Importance using the change in AUC\nwhen each predictor is permuted",
         cex.lab=3, cex.main=3, cex.axis=2)
    grid()  
    boxplot(t(xright[,-c(1:2)]), horizontal=TRUE, add=TRUE, at=ymiddle, yaxt="n", col="lightblue", xaxt="n")
    points(y=ymiddle, x=xright[,"train"], cex=3, pch=8, lwd=3, col="darkslateblue")
    points(y=ymiddle, x=xright[,"test"], cex=2.5, pch=5, lwd=3, col="royalblue")
    legend(x="bottomright", legend=c("CV","Train", "Test"),
           pch=c(22,8,5), pt.cex=c(3.5,2.75,2.5), pt.lwd=c(2,3,3),
           pt.bg=c("lightblue","darkslateblue", "royalblue"), col=c("black","darkslateblue", "royalblue"), cex=2.5)
  }
  
  Offset <- ifelse(!out$validationOptions$SplitData & !out$validationOptions$CrossValidate, 0.25, 0)
  ylabs <- rownames(xright)
  ylabs <- paste(substr(ylabs,start=1,stop=10), c("\n","")[1+(nchar(ylabs)<=10)], 
                 substr(ylabs,start=11,stop=nchar(ylabs)), sep="")
  axis(2, at=seq(from=0,to=length(out$inputVars), length=length(out$inputVars))+Offset, 
       labels=ylabs, las=2, cex=2.5, cex.lab=2.5, cex.axis=2.5)
  title(ylab="Variables",line=14, cex.lab=3, font.lab=2)
} 

### Confusion Matrix function --------------------------------------------------

confusion.matrix <- function(Stats,    # output from calcStat function 
                             out)      # out list
  {
  
  par(oma=c(4,3,5,3),mar=c(20,6,5,2))
  
  if(!out$validationOptions$SplitData & !out$validationOptions$CrossValidate){
    lo <- layout(matrix(data=c(1,2), nrow=1, ncol=2), c(4.5,1), 1)
  } else {
    lo <- layout(matrix(data=c(1,2,3), nrow=1, ncol=3), c(4.5,4.5,1), 1)

    if(out$validationOptions$CrossValidate){
      a <- lapply(Stats[-which(names(Stats) %in% c("train", "test"))],function(lst){lst$Cmx})
      cmx <- a[[1]]
      for(i in 2:length(a)){cmx <- cmx + a[[i]]} 
      csv.stats <- apply(do.call("rbind",(lapply(Stats[-which(names(Stats) %in% c("train", "test"))],
                                                 function(lst){
                                                   return(c(lst$Sens,lst$Specf,lst$Pcc,lst$Kappa,lst$Tss))
                                                   }))), 2, mean)
      Stats$crossValidation <- list(Cmx=cmx, Sens=csv.stats[1], Specf=csv.stats[2],
                                    Pcc=csv.stats[3], Kappa=csv.stats[4], Tss=csv.stats[5])
      Stats <- Stats[c("crossValidation", "train")]
      }
  }
  
  # zlim<-c(min(unlist(lapply(Stats,function(lst){100*lst$Cmx/sum(lst$Cmx)}))),max(unlist(lapply(Stats,function(lst){100*lst$Cmx/sum(lst$Cmx)}))))
  # instead of basing the zlim on the actual confusion matrices, base them on the maximum achievable value for a cell given the ratio of pres/abs
  
  extract.max <- function(lst){
    max(100*table(lst$auc.data$pres.abs)/length(lst$auc.dat$pres.abs))
  }
  options(warn=-1)
  zlim=c(0,100)
  options(warn=0)
  
  for(i in length(Stats):1){
    image((1:2),c(2,4),matrix(data=c(100*Stats[[i]]$Cmx[2]/sum(Stats[[i]]$Cmx[1:2]),
                                     100*Stats[[i]]$Cmx[4]/sum(Stats[[i]]$Cmx[3:4]),
                                     100*Stats[[i]]$Cmx[1]/sum(Stats[[i]]$Cmx[1:2]),
                                     100*Stats[[i]]$Cmx[3]/sum(Stats[[i]]$Cmx[3:4])),nrow=2),
          zlim=zlim,xaxt="n",yaxt="n",xlab="",
          ylab="",main=paste("Confusion matrix for \n", names(Stats)[i], "data",sep=" "),col=heat.colors(100)[100:1],cex.lab=2,cex.main=2.5)
    mtext("Absence",side=2,at=2,cex=2,lwd=1.3)
    mtext("Presence",side=2,at=4,cex=2,lwd=1.3)
    mtext("Presence",side=1,at=1,cex=2,line=1,lwd=1.3)
    mtext("Absence",side=1,at=2,cex=2,line=1,lwd=1.3)
    text(x=c(1,1,2,2),y=c(2,4,2,4),
         labels=c(Stats[[i]]$Cmx[2],Stats[[i]]$Cmx[1],
                  Stats[[i]]$Cmx[4],
                  Stats[[i]]$Cmx[3]),cex=5)
    abline(h=3,lwd=5)
    abline(v=1.5,lwd=5)
    mtext(paste(
      "Pct Correctly Classified : ",signif(Stats[[i]]$Pcc,digits=3),
      "\nSensitivity                      : ",signif(Stats[[i]]$Sens,digits=3),
      "\nSpecificity                      : ",signif(Stats[[i]]$Specf,digits=3),
      "\nTrue Skills Stat              : ",signif(Stats[[i]]$Tss,digits=3),
      "\nCohen's Kappa              : ",signif(Stats[[i]]$Kappa,digits=3),
      sep=""),
      side=1,line=13,cex=1.4,adj=0)
    box()
  }
  mtext("Observed",1,outer=TRUE,lwd=2,cex=2.5)
  mtext("Predicted",2,outer=TRUE,lwd=2,cex=2.5)
  
  ### color scale
  image(1,seq(from=zlim[1],to=zlim[2],length=50),
        matrix(data=seq(from=zlim[1],to=zlim[2],length=50), ncol=50,nrow=1),
        col=heat.colors(50)[50:1],
        xlab="",ylab="",zlim=zlim,
        xaxt="n",cex.lab=2,cex.axis=2)
  
}

### Test/Train ROC Plot function -----------------------------------------------

TestTrainRocPlot <- function(dat,    # Stats$train$auc.data
                             threshold = 101, find.auc = TRUE,
                             which.model = (1:(ncol(dat)-2)), na.rm = FALSE, 
                             xlab = "1-Specificity (false positives)",
                             ylab = "Sensitivity (true positives)", main = "ROC Plot",
                             model.names = NULL, color = NULL, line.type = NULL, lwd = 1,
                             mark = 0, mark.numbers = TRUE, mark.color = NULL, opt.thresholds = NULL,
                             opt.methods = NULL, req.sens, req.spec, obs.prev = NULL,
                             smoothing = 1, add.legend = TRUE, legend.text = model.names,
                             legend.cex = 0.8, add.opt.legend = TRUE, opt.legend.text = NULL,
                             opt.legend.cex = 0.7, counter.diagonal = FALSE, pch = NULL,
                             FPC, FNC, cost.line = FALSE,add.roc=FALSE)
{
  if (is.data.frame(dat) == FALSE) {
    if (is.matrix(dat) == TRUE) {dat <- as.data.frame(dat)
    } else {stop("'dat' must be either data frame or matrix")}}
  
  obs <- dat[, 2]
  if (length(obs[obs == 0]) == 0) { 
    stop("no observed absences in dataset, therefore specificity does not",
         "exist, and modeling, much less Area Under the Curve, is not very",
         "meaningful")}
  
  if (length(obs[obs == 1]) == 0) {
    stop("no observed presences in dataset, therefore sensitivity does not",
         "exist, and modeling, much less Area Under the Curve, is not very",
         "meaningful")}
  
  if (is.logical(find.auc) == FALSE) { 
    stop("'find.auc' must be of logical type")}
  
  if (is.logical(na.rm) == FALSE) {
    stop("'na.rm' must be of logical type")}
  
  if (is.logical(mark.numbers) == FALSE) {
    stop("'mark.numbers' must be of logical type!")}
  
  if (is.logical(add.legend) == FALSE) {
    stop("'add.legend' must be of logical type!")}
  
  if (is.logical(add.opt.legend) == FALSE) {
    stop("'add.opt.legend' must be of logical type")}
  
  if (is.logical(counter.diagonal) == FALSE) {
    stop("'counter.diagonal' must be of logical type")}
  
  if (length(smoothing) != 1) {
    stop("'smoothing' must be a single number greater than or equal to 1")
    } else {
    if (is.numeric(smoothing) == FALSE) {
      stop("'smoothing' must be a single number greater than or equal to 1")
    } else {
      if (smoothing < 1) {
        stop("'smoothing' must be a single number greater than or equal to 1")
      }}}
  
  if (sum(is.na(dat)) > 0) {
    if (na.rm == TRUE) {
      NA.rows <- apply(is.na(dat), 1, sum)
      warning(length(NA.rows[NA.rows > 0]), " rows ignored due to NA values")
      dat <- dat[NA.rows == 0, ]
    } else { return(NA) }}
  
  dat[dat[, 2] > 0, 2] <- 1
  N.models <- ncol(dat) - 2
  if (is.null(obs.prev) == TRUE) {
    obs.prev <- sum(dat[, 2])/nrow(dat)
    }
  
  if (obs.prev < 0 || obs.prev > 1) {
    stop("'obs.prev' must be a number between zero and one")} 
  
  if (obs.prev == 0) {
    warning("because your observed prevalence was zero, results may be strange")}
  
  if (obs.prev == 1) {
    warning("because your observed prevalence was one, results may be strange")}
  
  if (min(which.model) < 1 || sum(round(which.model) != which.model) != 0){
    stop("values in 'which.model' must be positive integers")}
  
  if (max(which.model) > N.models) {
    stop("values in 'which.model' must not be greater than number of models in 'dat'!")}
  
  if (is.null(model.names) == TRUE) {
    model.names <- if (is.null(names(dat)) == FALSE) {
      names(dat)[-c(1, 2)]
    } else {
      paste("Model", 1:N.models)
    }}
  
  if (N.models != length(model.names) && (length(which.model) != 1 || length(model.names) != 1)) {
    stop("If 'model.names' is specified it must either be a single name, or a vector",
         "of the same length as the number of model predictions in 'dat'")}
  
  if (is.null(legend.text) == TRUE) { legend.text <- model.names }
  
  if (length(legend.text) != N.models) {
    stop("'opt.legend.text' must be of same length as 'opt.methods'")}
  
  dat <- dat[, c(1, 2, which.model + 2)]
  if (length(model.names) != 1) { model.names <- model.names[which.model] }
  if (length(legend.text) != 1) { legend.text <- legend.text[which.model] }
  
  N.dat <- ncol(dat) - 2
  if (is.null(obs.prev) == TRUE) { obs.prev <- sum(dat[, 2])/nrow(dat) }
  if (obs.prev < 0 || obs.prev > 1) {
    stop("'obs.prev' must be a number between zero and one")}
  
  mark <- matrix(mark, length(mark), N.dat)
  if (!is.null(opt.methods) && is.null(opt.thresholds)) { opt.thresholds <- TRUE }
  if (is.null(opt.methods) && is.null(opt.thresholds)) { opt.thresholds <- FALSE }
  
  if (is.null(opt.methods)) { opt.methods <- c(1, 2, 4) }
  
  if (is.logical(opt.thresholds) == TRUE) {
    if (opt.thresholds == TRUE) {
      POSSIBLE.meth <- c("Default", "Sens=Spec", "MaxSens+Spec",
                         "MaxKappa", "MaxPCC", "PredPrev=Obs", "ObsPrev",
                         "MeanProb", "MinROCdist", "ReqSens", "ReqSpec",
                         "Cost")
      N.meth <- length(opt.methods)
      if (is.numeric(opt.methods) == TRUE) {
        if (sum(opt.methods %in% (1:length(POSSIBLE.meth))) != N.meth) {
          stop("invalid optimization method")
          } else {
          opt.methods <- POSSIBLE.meth[opt.methods]
          }}
      
      if (sum(opt.methods %in% POSSIBLE.meth) != N.meth) {
        stop("invalid optimization method") }
      
      if (is.null(opt.legend.text) == TRUE) { opt.legend.text <- opt.methods }
      
      if (length(opt.legend.text) != N.meth) {
        stop("'opt.legend.text' must be of same length as 'opt.methods'") }
      
      if ("ReqSens" %in% opt.methods) {
        if (missing(req.sens)) {
          warning("req.sens defaults to 0.85")
          req.sens <- 0.85
        }}
      
      if ("ReqSpec" %in% opt.methods) {
        if (missing(req.spec)) {
          warning("req.spec defaults to 0.85")
          req.spec <- 0.85
        }}
      
      if ("Cost" %in% opt.methods) {
        if (missing(FPC) || missing(FNC)) {
          warning("costs assumed to be equal")
          FPC <- 1
          FNC <- 1
        }
        if (FPC <= 0 || FNC <= 0) { stop("costs must be positive") }
        if (is.logical(cost.line) == FALSE) {
          stop("'cost.line' must be of logical type")
          }
        if (!"Cost" %in% opt.methods) { cost.line <- FALSE }
      }
      
      mark <- optimal.thresholds(dat = dat, threshold = threshold,
                                 model.names = model.names, na.rm = na.rm, 
                                 opt.methods = opt.methods, req.sens = req.sens, 
                                 req.spec = req.spec, obs.prev = obs.prev,
                                 smoothing = smoothing, FPC = FPC, FNC = FNC)[,-1, drop = FALSE]
      
      if (is.null(pch) == TRUE) {
        pch <- c(1, 5, 2, 16, 15, 17, 8, 6, 9, 12, 4,
                 7)[match(opt.methods, POSSIBLE.meth)]
      } else {
        pch <- rep(pch, length(opt.methods))[1:length(opt.methods)]
      }
    }
  }
  
  if (is.logical(opt.thresholds) == FALSE) {
    if (!is.numeric(opt.thresholds)) {
      stop("'opt.thresholds' must be 'TRUE', 'FALSE', or numeric") }
    
    if (min(opt.thresholds) < 0) { 
      stop("'opt.thresholds' can not be negative") }
    
    if (max(opt.thresholds) > 1) {
      if (N.thr == 1 && round(opt.thresholds) == opt.thresholds) {
        opt.thresholds <- seq(length = opt.thresholds,
                              from = 0, to = 1)
        N.thr <- length(opt.thresholds)
      } else {
        stop("non-interger, non-logical 'opt.thresholds' greater than 1")
      }}
    
    N.opt.thresh <- length(opt.thresholds)
    if (is.null(opt.legend.text)) {
      opt.legend.text <- rep("threshold", N.opt.thresh)
    }
    
    if (length(opt.legend.text) != N.opt.thresh) {
      stop("length of 'opt.legend.text' does not match number of specified thresholds") }
    
    if (is.null(pch)) { pch <- 1:N.opt.thresh }
    
    if (length(pch) != N.opt.thresh) {
      stop("length of 'pch' does not match number of specified thresholds") }
    
    mark <- matrix(opt.thresholds, length(opt.thresholds), N.dat)
    opt.thresholds = TRUE
  }
  
  if (is.null(pch) == TRUE) { pch <- 16 }
  
  if (is.null(color) == TRUE) { 
    colors <- rep(1, N.dat)
    if (is.null(line.type) == TRUE) {
      line.type <- (1:N.dat) + 1
    }} else {
    if (is.logical(color) == TRUE) {
      if (color == FALSE) {
        colors <- rep(1, N.dat)
        if (is.null(line.type) == TRUE) {
          line.type <- (1:N.dat) + 1
        }} else {
        colors <- (1:N.dat) + 1
        lwd <- 2 * lwd
        if (is.null(line.type) == TRUE) {
          line.type <- rep(1, N.dat)
        }}} else {
          colors <- rep(color, N.dat)[1:N.dat]
          lwd <- 2 * lwd
          if (is.null(line.type) == TRUE) {
            line.type <- rep(1, N.dat)
          }}}
  if (is.null(mark.color) == TRUE) {
    mark.colors <- colors
  } else {
    if (is.logical(mark.color) == TRUE) {
      if (mark.color == FALSE) {
        mark.colors <- rep(1, N.dat)
      } else {
        mark.colors <- (1:N.dat) + 1
      }} else {
      mark.colors <- rep(mark.color, N.dat)[1:N.dat]
      }}
  if (is.null(line.type) == FALSE) {
    if (is.logical(line.type) == TRUE) {
      if (line.type == FALSE) {
        line.type <- rep(1, N.dat)
      } else {
        line.type <- (1:N.dat) + 1
      }} else {
      line.type <- rep(line.type, N.dat)[1:N.dat]
    }}
  op <- par(pty = "s")
  if(add.roc==FALSE){
    plot(c(0, 1), c(0, 1), type = "n", xlab = xlab, ylab = ylab,
         main = main,cex.main=2.2,cex.lab=2,cex.axis=1.5)}
  lines(c(0, 1), c(0, 1), col = "lightgray")
  if (counter.diagonal == TRUE) {
    abline(a = 1, b = -1, col = "lightgray")
  }
  for (d in 1:N.dat) {
    Model.dat <- roc.plot.calculate(DATA = dat, 
                                    threshold = threshold,
                                    which.model = d)
    lines(x = (1 - Model.dat$specificity), y = Model.dat$sensitivity,
          lty = line.type[d], lwd = lwd, col = colors[d])
    if (max(mark) != 0) {
      Mark.dat <- roc.plot.calculate(DATA = dat, 
                                     threshold = mark[,d], 
                                     which.model = d)
      
      Mark.pretty <- round(Mark.dat$threshold, 2)
      Mark.pretty.char <- as.character(Mark.pretty)
      Mark.pretty.char[Mark.pretty == 0] <- "0.00"
      Mark.pretty.char[Mark.pretty == 1] <- "1.00"
      Mark.pretty.char[nchar(Mark.pretty.char) == 3] <- paste(Mark.pretty.char[nchar(Mark.pretty.char) ==
                                                                                 3], "0", sep = "")
      if(add.roc==FALSE){
        points(x = (1 - Mark.dat$specificity), y = Mark.dat$sensitivity,
               cex = 2, pch = pch, col = mark.colors[d])
        
        if (mark.numbers == TRUE) {
          text(x = (1 - Mark.dat$specificity), y = Mark.dat$sensitivity -
                 0.03, labels = Mark.pretty.char, pos = 4, col = mark.colors[d])
        }}
    }
    if (cost.line == TRUE) {
      tag <- match("Cost", opt.methods)
      sl <- (FPC/FNC) * (1 - obs.prev)/obs.prev
      if (obs.prev == 0) {
        obs.prev <- 1e-06
      }
      abline(a = Mark.dat$sensitivity[tag] - ((1 - Mark.dat$specificity[tag]) *
                                                sl), b = sl, col = mark.colors[d], lty = 3)
    }
  }
  
  inset <- c(0.02, 0.02)
  if (opt.thresholds == TRUE && add.opt.legend == TRUE) {
    if (N.dat == 1) {
      opt.legend.names <- paste(Mark.pretty.char, opt.legend.text)
    }
    else {
      opt.legend.names <- opt.legend.text
    }
    leg.loc <- legend(x = "bottomright", inset = inset, pt.cex = 1,
                      legend = opt.legend.names, pch = pch, bg = "white",
                      cex = 1.3)
    inset <- c(0.02, leg.loc$rect$top + 0.05)
  }
  if (add.legend == TRUE) {
    legend.names <- legend.text
    if (find.auc == TRUE) {
      AUC <- rep(0, N.dat)
      for (dat in 1:N.dat) {
        AUC[dat] <- auc(dat = dat, which.model = dat)$AUC
      }
      AUC.pretty <- round(AUC, 2)
      AUC.pretty.char <- as.character(AUC.pretty)
      AUC.pretty.char[AUC.pretty == 0] <- "0.00"
      AUC.pretty.char[AUC.pretty == 1] <- "1.00"
      AUC.pretty.char[nchar(AUC.pretty.char) == 3] <- paste(AUC.pretty.char[nchar(AUC.pretty.char) ==
                                                                              3], "0", sep = "")
      legend.names <- paste(AUC.pretty.char, legend.text)
    }
    legend(x = "bottomright", inset = inset, legend = legend.names,
           lty = line.type, col = colors, title = "AUC:", cex = legend.cex,
           lwd = lwd, bg = "white")
  }
  par(op)
}

### Presence-Only Smoothed Calibration Plot function ---------------------------

pocplot <- function(pred, back, linearize=TRUE, ...){
  
  ispresence <- c(rep(1,length(pred)), rep(0, length(back)))
  predd <- smoothdist(c(pred,back), ispresence)
  c <- mean(back)*length(back)/length(pred)
  if (linearize) {
    fun <- function(x,y) c*y / (1-y)
    predd$y <- mapply(fun, predd$x, predd$y)
    predd$se <- mapply(fun, predd$x, predd$se)
    ideal <- function(x) x
    ylab <- "Relative probability of presence" 
  } 
  else {
    ideal <- function(x) x / (x + c)
    ylab <- "Probability of presence"
  }
  calibplot(predd, negrug=back, posrug=pred, ideal=ideal, ylab=ylab,
            capuci = FALSE, ...)
  predd
}

### Presence-Absence Smoothed Calibration Plot function ------------------------

pacplot <- function(pred, pa, ...) {
  predd <- smoothdist(preds = pred, obs = pa)
  calibplot(predd, negrug=pred[pa==0], posrug=pred[pa==1], 
            ideal=function(x) x, ylab="Probability of presence", ...)
}

#### Plotting function for Calibration plots [nested in pocplot/pacplot] -------

calibplot <- function(pred, negrug, posrug, ideal, 
                      ylim=c(0,1), capuci=TRUE,
                      xlabel = "Predicted probability of presence", 
                      filename=NULL, title="Calibration plot", ...){
  
  if (!is.null(filename)){ png(filename) }
  ylow <- pred$y - 2 * pred$se
  ylow[ylow<0] <- 0
  yhigh <- pred$y + 2 * pred$se
  if (capuci) yhigh[yhigh>1] <- 1
  plot(pred$x, ylow, type="l", col="grey", 
       ylim=ylim, xlim=range(pred$x), main=title,
       xlab=xlabel, lwd=2,cex.axis=1.4,cex.lab=1.8,cex.main=2, ...)
  lines(pred$x, yhigh, lwd=2, col="grey")
  lines(pred$x, sapply(pred$x, ideal), lty="dashed")
  points(pred$x, pred$y, col="blue")
  rug(negrug,col="blue",lwd=2)
  rug(posrug, col = "red",lwd=2)
  
  if (!is.null(filename)){ dev.off() }
}

#### Smoothing function for Calibration plots [nested in pocplot/pacplot] ------

smoothingdf <- 6
smoothdist <- function(preds, obs) {
  
  gam1 <- try(glm(obs ~ ns(preds, df=smoothingdf), 
                  weights=rep(1, length(preds)), 
                  family=binomial),silent=TRUE)
  
  if("try-error" %in% class(gam1)){
    gam1 <- try(glm(obs ~ ns(preds, df=1), 
                    weights=rep(1, length(preds)), 
                    family=binomial),silent=TRUE)
  }
  
  x <- seq(min(preds), max(preds), length = 512)
  y <- predict(gam1, newdata = data.frame(preds = x), se.fit = TRUE,
               type = "response")
  data.frame(x=x, y=y$fit, se=y$se.fit)
}

### Capture Statistics function ------------------------------------------------


capture.stats <- function(Stats.lst,  # stats or lst output from calcStat function 
                          file.name,
                          label,      # splitType
                          family,
                          opt.methods,
                          out){
  
  if(label == "train/test"){ label = "Final evaluation" }
  
  capture.output(cat(" applied to",label, "split:\n",sep=" "),
                 file=file.name,append=TRUE)
  capture.output(cat( "\n",
                      "\n\t Correlation Coefficient      :",mean(unlist(lapply(Stats.lst,function(lst){lst$correlation}))),
                      if(label %in% c("cv", "cv/train/test")){
                        paste(" (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$correlation}))),digits=5), ")",sep="")},
                      if(!out$hasSplit |(out$hasSplit & label == "none")){ 
                        paste("\n\t NULL Deviance                : ", signif(mean(unlist(lapply(Stats.lst,function(lst){lst$null.dev}))),digits=5),
                        if(out$hasSplit & label == "none"){" (Averaged over background splits)"},
                        if(label %in% c("cv", "cv/train/test")){
                          paste(" (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$null.dev}))),digits=5), ")",sep="")},
                        "\n\t Fit Deviance                 : ", signif(mean(unlist(lapply(Stats.lst,function(lst){lst$dev.fit}))),digits=5),
                        if(out$hasSplit & label=="none"){" (Averaged over background splits)"},  
                        if(label %in% c("cv", "cv/train/test")){
                          paste(" (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$dev.fit}))),digits=5),")",sep="")},
                        "\n\t Explained Deviance           : ", signif(mean(unlist(lapply(Stats.lst,function(lst){lst$dev.exp}))),digits=5),
                        if(label %in% c("cv", "cv/train/test")){
                          paste(" (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$dev.exp}))),digits=5),")",sep="")},
                        "\n\t Percent Deviance Explained   : ", signif(mean(unlist(lapply(Stats.lst,function(lst){lst$pct.dev.exp}))),digits=5),
                        if(label %in% c("cv", "cv/train/test")){
                          paste(" (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$pct.dev.exp}))),5), ")",sep="")},sep="")},
                      file=file.name, append=TRUE))
  
  if(family %in% c("binomial","bernoulli")){
    capture.output(cat(
      "\n\n  Threshold Methods based on", switch(opt.methods,
                                                 "1"=".5 threshold",
                                                 "2"="Sens=Spec",
                                                 "3"="maximize (sensitivity+specificity)/2",
                                                 "4"="maximize Kappa",
                                                 "5"="maximize percent correctly classified",
                                                 "6"="predicted prevalence=observed prevalence",
                                                 "7"="threshold=observed prevalence",
                                                 "8"="mean predicted probability",
                                                 "9"="minimize distance between ROC plot and (0,1)"),
      if (label %in% c("none", "Final evaluation")){
        paste("\n\t Threshold                    : ", Stats.lst[[1]]$thresh)
        } else { 
          paste("\n\t Mean Threshold               : ", mean(unlist(lapply(Stats.lst,function(lst){lst$thresh}))),
        " (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$thresh}))),digits=5),")",sep="")
      },
      "\n\n\t Confusion Matrix: \n\n"),
      if(label %in% c("none", "Final evaluation")){ print.table(Stats.lst[[1]]$Cmx)
      } else {
        a <-lapply(Stats.lst,function(lst){lst$Cmx})
        cmx <-a[[1]]
        for(i in 2:length(a)) cmx<-cmx+a[[i]] #it's amazing I can't think of a better way to sum a list of tables
        print.table(cmx)
      },
      cat(
        "\n\t AUC                          : ",mean(unlist(lapply(Stats.lst,function(lst){lst$auc.fit}))),
        if(label %in% c("cv", "cv/train/test")){
          paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$auc.fit}))),digits=5),")",sep="")},
        "\n\t AUC-pr                       : ",mean(unlist(lapply(Stats.lst,function(lst){lst$auc.pr}))),
        if(label %in% c("cv", "cv/train/test")){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$auc.pr}))),digits=5),")",sep="")},
        
        "\n\t Percent Correctly Classified : ",mean(unlist(lapply(Stats.lst,function(lst){lst$Pcc}))),
        if(label %in% c("cv", "cv/train/test")){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$Pcc}))),digits=5),")",sep="")},
        
        "\n\t Sensitivity                  : ",mean(unlist(lapply(Stats.lst,function(lst){lst$Sens}))),
        if(label %in% c("cv", "cv/train/test")){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$Sens}))),digits=5),")",sep="")},
        
        "\n\t Specificity                  : ",mean(unlist(lapply(Stats.lst,function(lst){lst$Specf}))),
        if(label %in% c("cv", "cv/train/test")){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$Specf}))),digits=5),")",sep="")},
        
        "\n\t Kappa                        : ",mean(unlist(lapply(Stats.lst,function(lst){lst$Kappa}))),
        if(label %in% c("cv", "cv/train/test")){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$Kappa}))),digits=5),")",sep="")},
        
        "\n\t True Skill Statistic         : ",mean(unlist(lapply(Stats.lst,function(lst){lst$Tss}))),
        if(label %in% c("cv", "cv/train/test")){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$Tss}))),digits=5),")",sep="")},
        "\n"),
      
      file=paste0(out$tempDir, out$modType,"_output.txt",sep=""),append=TRUE)
  }
  
  if(!out$pseudoAbs){
    capture.output(cat( "\n\n   Calibration Statistics",
                                             "\n\t Intercept (general calibration)                            : ",mean(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[1]}))),
                                             if(label %in% c("cv", "cv/train/test")){paste(" (sd ",
                                                                                signif(sd(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[1]}))),digits=5), ")",sep="")},
                                             "\n\t Slope   (direction and variation in fit)                   : ",mean(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[2]}))),
                                             if(label %in% c("cv", "cv/train/test")){paste(" (sd ",
                                                                                signif(sd(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[2]}))),digits=5),
                                                                                ")",sep="")},
                                             "\n\t Testa0b1 (overall reliability of predictors)               : ",mean(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[3]}))),
                                             if(label %in% c("cv", "cv/train/test")){paste(" (sd ",
                                                                                signif(sd(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[3]}))),digits=5),
                                                                                ")",sep="")},
                                             "\n\t Testa0|b1(incorrect calibration given correct refinement)  : ",mean(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[4]}))),
                                             if(label %in% c("cv", "cv/train/test")){paste(" (sd ",
                                                                                signif(sd(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[4]}))),digits=5),
                                                                                ")",sep="")},
                                             "\n\t Testb1|a (refinement given correct calibration)            : ",mean(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[5]}))),
                                             if(label %in% c("cv", "cv/train/test")){paste(" (sd ",
                                                                                signif(sd(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[5]}))),digits=5),
                                                                                ")",sep="")},
                                             
                                             "\n\n",
                                             file=paste0(out$tempDir, out$modType, "_output.txt",sep=""),append=TRUE))
  #if(label=="crossValidation"){cat("\n\n   Pooled Calibration Statistics\n",print.table(cbind(names(out$cv$pooled.calib),out$cv$pooled.calib)))}
  #something I should include later
  }
}
