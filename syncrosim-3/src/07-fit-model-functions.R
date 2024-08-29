## ---------------------------- 
## wisdm - fit model functions  
## ApexRMS, March 2022       
## ---------------------------- 

library(PresenceAbsence)
library(PRROC)
library(ROCR)
library(ggplot2)
library(splines)

# MODEL FIT FUNCTION -----------------------------------------------------------

## Fit model -------------------------------------------------------------------

fitModel <- function(dat,           # df of training data 
                     out,           # out list
                     fullFit=TRUE,
                     pts=NULL,
                     weight=NULL){
                     # Fold,...
                     
  # This function was written to separate the steps involved in model fitting 
  # from the post processing steps needed to produce several of the later outputs
  # so the same function call can be used for cross-validation and generic model fit.  
  
  out$seed <- 32639
  set.seed(out$seed)
  
  # Sanitize output variable names
  sanitizedVarNames <- out$inputVars
  for (i in 1:length(sanitizedVarNames)) {
    if (grepl("\\s", sanitizedVarNames[i])){
      oldVarName <- sanitizedVarNames[i]
      newVarName <- paste0("`", oldVarName, "`")
      sanitizedVarNames[i] <- newVarName
    }
  }
  
  #================================================================
  #                        GLM
  #================================================================= 
  
  if(out$modType == "glm") {
    
    if(out$pseudoAbs){
      # for glm sum of absence weights set equal to sum of presence weights
      absWt <- sum(dat$Response == 1)/sum(dat$Response == 0)
      dat$Weight[dat$Response == 0] <- absWt
    }
    
    penalty <- if(out$modOptions$SimplificationMethod == "AIC"){2} else {log(nrow(dat))}
    
    factor.mask <- na.omit(match(out$factorInputVars,sanitizedVarNames))
    cont.mask <- seq(1:length(out$inputVars))
    if(length(factor.mask)!=0){cont.mask<-cont.mask[-c(factor.mask)]}
    
    if(!out$modOptions$ConsiderSquaredTerms & !out$modOptions$ConsiderInteractions){
      scopeGLM <- list(lower = as.formula(paste("Response","~1")),
                       upper = as.formula(paste("Response","~",paste(sanitizedVarNames, collapse='+'))))
    }
    if(out$modOptions$ConsiderSquaredTerms & out$modOptions$ConsiderInteractions){ # creates full scope with interactions and squared terms
      scopeGLM <- list(lower = as.formula(paste("Response","~1")),
                       upper = as.formula(paste("Response","~",paste(c(if(length(factor.mask)>0) paste(sanitizedVarNames[factor.mask],collapse=" + "),
                                                                       paste("(",paste(sanitizedVarNamess[cont.mask],collapse=" + "),")^2",sep=""),
                                                                       paste("I(",sanitizedVarNames[cont.mask],"^2)",sep="")),collapse=" + "),sep="")))
    }
    if(!out$modOptions$ConsiderSquaredTerms & out$modOptions$ConsiderInteractions){ # creates full scope with interactions
      scopeGLM <- list(lower = as.formula(paste("Response","~1")),
                       upper = as.formula(paste("Response","~",paste(c(if(length(factor.mask)>0) paste(sanitizedVarNames[factor.mask],collapse=" + "),
                                                                       paste("(",paste(sanitizedVarNames[cont.mask],collapse=" + "),")^2",sep="")),
                                                                     collapse=" + "),sep="")))
    }
    if(out$modOptions$ConsiderSquaredTerms & !out$modOptions$ConsiderInteractions){ # creates full scope with squared terms
      scopeGLM <- list(lower = as.formula(paste("Response","~1")),
                       upper = as.formula(paste("Response","~",paste(c(paste(sanitizedVarNames,collapse=" + "),
                                                                       paste("I(",sanitizedVarNames[cont.mask],"^2)",sep="")),collapse=" + "),sep="")))
    }
    
    if(out$modOptions$SelectBestPredictors){
      modelGLMStep <- step(glm(scopeGLM$lower,family=out$modelFamily,data=dat,weights=dat$Weight,na.action="na.exclude"),
                             direction='both',scope=scopeGLM,k=penalty,trace=1)
    } else {
      modelGLMStep <- glm(scopeGLM$upper,family=out$modelFamily,data=dat,weights=dat$Weight,na.action="na.exclude")
    }
    return(modelGLMStep)
  }


  #================================================================
  #                        RF
  #=================================================================
  if(out$modType == "rf"){
    
    # set defaults
    n.trees = out$modOptions$NumberOfTrees
    nodesize = out$modOptions$NodeSize
    importance = out$modOptions$EvaluateCovariateImportance
    localImp = out$modOptions$CalculateCasewiseImportance
    samp.replace = out$modOptions$SampleWithReplacement
    norm.votes = out$modOptions$NormalizeVotes
    proximity = out$modOptions$CalculateProximity 
    
    if(is.na(out$modOptions$NumberOfVariablesSampled)){
      mtry = NULL } else { mtry = out$modOptions$NumberOfVariablesSampled }
    if(is.na(out$modOptions$MaximumNodes)){
      maxnodes = NULL } else { maxnodes = out$modOptions$MaximumNodes }

    sampsize=NULL
    nPerm=1
    oob.prox=proximity
    keep.forest=NULL
    keep.inbag=FALSE
    xtest=NULL
    ytest=NULL

    if("predicted" %in% names(dat)){ dat <- select(dat, -predicted)}
    
    x = dat[,-c(1:7)] 
    if(out$modelFamily == "poisson"){ 
      y = dat$Response
    } else { y = factor(dat$Response) }
      
    # tune the mtry parameter - this controls the number of covariates randomly subset for each split #
    if(is.null(mtry)){
      cat("\ntuning mtry parameter\n")
      mtr <- tuneRF(x = x, y = y,
                    mtryStart = 3,
                    importance = TRUE,
                    ntreeTry = 100,
                    replace=FALSE, 
                    doBest=F,
                    plot=F)
      mtry <- mtr[mtr[,2] == min(mtr[,2]),1][1]
    }
    
    cat("\nnow fitting full random forest model using mtry=",mtry,"\n")
    modelRF <- randomForest(
      x = x,
      y = y,
      xtest = xtest,
      ytest = ytest,
      importance = importance,
      ntree = n.trees,
      mtry = mtry,
      replace = samp.replace,
      sampsize = ifelse(is.null(sampsize),
                        rep(min(ceiling(table(y) * 0.632)), 2), # balanced downsample based on smallest class
                        sampsize),
      nodesize = nodesize,
      # ifelse(is.null(nodesize),(if (!is.null(y) && !is.factor(y)) 5 else 1),nodesize),
      maxnodes = maxnodes,
      localImp = localImp,
      nPerm = nPerm,
      keep.forest = ifelse(
        is.null(keep.forest),
        !is.null(y) && is.null(xtest),
        keep.forest
      ),
      corr.bias = corr.bias,
      keep.inbag = keep.inbag
    )
  return(modelRF)
  }


  #================================================================
  #                        MAXENT
  #=================================================================
  if(out$modType == "maxent"){
    
      # If there are parentheses in the working folder name, then the jar file will not run properly
      validTempPath <- sum(
        grepl("\\(", strsplit(out$tempDir, "")[[1]])
        ) + sum(
          grepl("\\)", strsplit(out$tempDir, "")[[1]])
          ) == 0
      
      if (!validTempPath){
        stop(paste0("Maxent model will not run if there are parentheses in the",
                    " temporary directory path. Please set a different ",
                    "temporary folder before continuing."))
      }

      if(fullFit == T){
      
      # prepare batch file
      capture.output(cat("java -mx", out$modOptions$MemoryLimit, "m", sep=""), file = out$batchPath)
      cat(" -jar", paste0('"', file.path(ssimEnvironment()$PackageDirectory, "maxent.jar", fsep = "\\"), '"'), file = out$batchPath, append = T)
      if(!out$modOptions$VisibleInterface){ cat(" -z", file = out$batchPath, append = T) }
      cat(paste0(' samplesfile="',out$swdPath, '"'), file = out$batchPath, append = T)
      cat(paste0(' environmentallayers="', out$backgroundPath, '"'), file = out$batchPath, append = T)
      if(length(out$factorInputVars)>0){
        for (i in out$factorInputVars){
          cat(paste0(" togglelayertype=",i), file = out$batchPath, append = T)
        }}
      if(!is.null(out$testDataPath)){
        cat(paste0(' testsamplesfile="',out$testDataPath,'"'), file = out$batchPath, append = T)
      }
      cat(paste0(' outputdirectory="',file.path(out$tempDir, "Outputs", fsep = "\\"),'"'), file = out$batchPath, append = T)
      cat(" threads=",out$modOptions$MultiprocessingThreads, sep = "", file = out$batchPath, append = T)
      cat(" responsecurves jackknife writeclampgrid writemess warnings prefixes", file = out$batchPath, append = T) # reverse these default settings  
      cat(" redoifexists autorun", file = out$batchPath, append = T)
      
      # Note than maxent can't handle spaces in the batch file path
      # - if there are spaces in tempDir, copy the batch file to a system temp file
      # - also update batchPath location in the local scope
      if (str_detect(out$tempDir, " ")) {
        batchTempFile <- tempfile(pattern = "runMaxent", fileext = ".bat")
        file.copy(out$batchPath, batchTempFile, overwrite = T)
        out$batchPath <- batchTempFile
      }
      # run maxent
      shell(out$batchPath)
      
      # read lambdas output
      modelMaxent <- read.maxent(file.path(out$tempDir, "Outputs", "species.lambdas", fsep = "\\"))
      
    } else {
      
      # prepare batch file
      capture.output(cat("java -mx", out$modOptions$MemoryLimit, "m", sep=""), file = out$batchPath)
      cat(" -jar", paste0('"', file.path(ssimEnvironment()$PackageDirectory, "maxent.jar", fsep = "\\"), '"'), file = out$batchPath, append = T)
      if(!out$modOptions$VisibleInterface){ cat(" -z", file = out$batchPath, append = T) }
      cat(paste0(' samplesfile="',file.path(out$tempDir, "CVsplits", "training-swd.csv", fsep = "\\"), '"'), file = out$batchPath, append = T)
      cat(paste0(' environmentallayers="', file.path(out$tempDir, "CVsplits", "background-swd.csv", fsep = "\\"), '"'), file = out$batchPath, append = T)
      if(length(out$factorInputVars)>0){
        for (i in out$factorInputVars){
          cat(paste0(" togglelayertype=",i), file = out$batchPath, append = T)
        }}
      cat(paste0(' outputdirectory="',file.path(out$tempDir, "CVsplits", fsep = "\\"),'"'), file = out$batchPath, append = T)
      cat(" threads=",out$modOptions$MultiprocessingThreads, sep = "", file = out$batchPath, append = T)
      cat(" writeclampgrid writemess warnings prefixes", file = out$batchPath, append = T) # reverse these default settings  
      cat(" redoifexists autorun", file = out$batchPath, append = T)
      
      # Note than maxent can't handle spaces in the batch file path
      # - if there are spaces in tempDir, copy the batch file to a system temp file
      # - also update batchPath location in the local scope
      if (str_detect(out$tempDir, " ")) {
        batchTempFile <- tempfile(pattern = "runMaxent", fileext = ".bat")
        file.copy(out$batchPath, batchTempFile, overwrite = T)
        out$batchPath <- batchTempFile
      }
      
      # run maxent
      shell(out$batchPath)
      
      # read lambdas output
      modelMaxent <- read.maxent(file.path(out$tempDir, "CVsplits", "species.lambdas", fsep = "\\"))
    
    }
    return(modelMaxent)
  }
  #================================================================
  #                        BRT
  #=================================================================
  if(out$modType == "brt"){
    
    # set n-folds
    if(out$validationOptions$CrossValidate){
      nFolds <- out$validationOptions$NumberOfFolds
    } else { nFolds <- 10 }
      
    # calculating the case weights
    prNum <- as.numeric(table(dat$Response)["1"])                                # number of presences
    bgNum <- as.numeric(table(dat$Response)["0"])                                # number of backgrounds
    wt <- ifelse(dat$Response == 1, 1, prNum / bgNum)
    
    # tmp <- Sys.time()
    if(fullFit){
      modelBRT <- gbm.step(data = dat,
                           gbm.x = 8:ncol(dat),                                  # column indices for covariates
                           gbm.y = 4,                                            # column index for response
                           family = out$modelFamily,
                           tree.complexity = ifelse(prNum < 50, 1, 5),
                           learning.rate = out$modOptions$LearningRate,
                           bag.fraction = out$modOptions$BagFraction,
                           max.trees = out$modOptions$MaximumTrees,
                           n.trees = out$modOptions$NumberOfTrees,
                           # step.size = out$modOptions$stepSize,
                           n.folds = nFolds,                                     # number of cross-validation folds
                           site.weights = wt,
                           plot.main = FALSE)                                    # avoid plotting hold-out deviance curve
                           # silent = TRUE)                                       # avoid printing the cv results
      
    } else {
      modelBRT <- gbm.step(data = dat,
                           gbm.x = 8:ncol(dat),                                  # column indices for covariates
                           gbm.y = 4,                                            # column index for response
                           family = out$modelFamily,
                           tree.complexity = ifelse(prNum < 50, 1, 5),
                           learning.rate = out$modOptions$LearningRate,
                           bag.fraction = out$modOptions$BagFraction,
                           max.trees = out$modOptions$MaximumTrees,
                           n.trees = out$modOptions$NumberOfTrees,
                           # step.size = out$modOptions$stepSize,
                           n.folds = nFolds,                                     # number of cross-validation folds
                           site.weights = wt,
                           plot.main = FALSE)                                    # avoid plotting hold-out deviance curve
      
    }

    
     # Sys.time() - tmp
     return(modelBRT)
  }
    # brt.full<-list()
    # lr.list<-list()
    # mod.simp<-list()
    # 
    # 
    # if(!is.null(tc)) out$mods$parms$tc.full<-out$mods$parms$tc.sub<-tc
    # 
    # #going to try to estimate learning rate and predictors to use in final model not just on the subset but by calculating for
    # #several of the splits (if the used was split)
    # lr.samp<-sample(1:num.splits,size=min(num.splits,5),replace=FALSE)
    # for(i in 1:length(lr.samp)){
    #   if(length(lr.samp)>1) {out$dat$Subset$dat<-dat[c(Split,rep(lr.samp[i],times=sum(dat$response>0)))==lr.samp[i],]
    #   out$dat$Subset$weight<-weight[c(Split,rep(lr.samp[i],times=sum(dat$response>0)))==lr.samp[i]]
    #   if(is.null(out$dat$Subset$weight))  out$dat$Subset$weight<-rep(1,times=nrow(out$dat$Subset$dat))
    #   out$dat$Subset$ratio=.5
    #   }
    #   lr.list[[i]]<-est.lr(out)
    # }
    # #now reassembling everything from lr estimation before continuing
    # out$mods$lr.mod$good.cols<-unique(unlist(lapply(lr.list,function(lst){lst$lr.mod$good.cols})))
    # out$mods$parms$tc.sub<-round(mean(unlist(lapply(lr.list,function(lst){lst$parms$tc.sub}))))
    # out$mods$parms$tc.full<-round(mean(unlist(lapply(lr.list,function(lst){lst$parms$tc.full}))))
    # out$mods$lr.mod$lr0<-mean(unlist(lapply(lr.list,function(lst){lst$lr.mod$lr0})))
    # out$mods$lr.mod$lr<-mean(unlist(lapply(lr.list,function(lst){lst$lr.mod$lr})))
    # 
    # cat("\nfinished with learning rate estimation, lr=",out$mods$lr.mod$lr0)
    # cat("\nfor final fit, lr=",out$mods$lr.mod$lr,"and tc=",out$mods$parms$tc.full,"\n")
    # 
    # if(simp.method=="cross-validation"){
    #   for(i in 1:length(lr.samp)){
    #     if(length(lr.samp)>1) {out$dat$Subset$dat<-dat[c(Split,rep(lr.samp[i],times=sum(dat$response>0)))==lr.samp[i],]
    #     out$dat$Subset$weight<-weight[c(Split,rep(lr.samp[i],times=sum(dat$response>0)))==lr.samp[i]]
    #     out$dat$Subset$ratio=.5
    #     }
    #     # remove variables with <1% relative influence and re-fit model
    #     if(length(out$mods$lr.mod$good.cols)<=1) stop("BRT must have at least two independent variables")
    #     max.trees<-NULL
    #     #learning rate estimation removes columns with low contributions to fit for removal
    #     #here we put specify use all if no predictor selection was to occur
    #     if(!predSelect) out$mods$lr.mod$good.cols<-seq(from=2,to=ncol(out$dat$Subset$dat))
    # 
    #     n.trees<-c(300,600,800,1000,1200,1500,1800)
    #     if(!is.null(out$input$n.trees)) n.trees=out$input$n.trees
    #     m0 <- gbm.step.fast(dat=out$dat$Subset$dat,gbm.x=out$mods$lr.mod$good.cols,gbm.y=1,family=model.family,
    #                         n.trees = n.trees,step.size=step.size,max.trees=max.trees,
    #                         tolerance.method=tolerance.method,tolerance=tolerance, n.folds=n.folds,prev.stratify=prev.stratify,
    #                         tree.complexity=out$mods$parms$tc.sub,learning.rate=out$mods$lr.mod$lr0,bag.fraction=bag.fraction,site.weights=out$dat$Subset$weight,
    #                         autostop=T,debug.mode=F,silent=!debug.mode,
    #                         plot.main=F,superfast=F)
    #     if(predSelect) mod.simp[[i]] <- gbm.simplify(m0,n.folds=n.folds,plot=F,verbose=F,alpha=alpha) # this step is very slow #
    # 
    #   }        #if we removed bad predictors the good predictor list otherwise make sure we specify include all again
    #   if(predSelect) out$mods$simp.mod$good.cols <- unique(unlist(lapply(mod.simp,function(lst){lst$pred.list[[length(lst$pred.list)]]})))
    #   else out$mods$simp.mod$good.cols <- seq(from=2,to=ncol(out$dat$Subset$dat))
    #   out$mods$simp.mod$good.vars <- names(dat)[out$mods$simp.mod$good.cols]
    #   {cat("\n");cat("50%\n")}
    # }
    # 
    # final.mod<-list()
    # 
    # for(i in 1:num.splits){
    #   if(out$mods$lr.mod$lr==0) out$mods$lr.mod$lr<-out$mods$lr.mod$lr0
    #   n.trees<-c(300,600,800,1000,1200,1500,1800,2200,2600,3000,3500,4000,4500,5000)
    #   if(!is.null(out$input$n.trees)) n.trees=out$input$n.trees
    #   final.mod[[i]] <- gbm.step.fast(dat=dat[c(Split,rep(i,times=sum(dat$response>0)))==i,],gbm.x=out$mods$simp.mod$good.cols,gbm.y = 1,family=model.family,
    #                                   n.trees = n.trees,n.folds=n.folds,max.trees,
    #                                   tree.complexity=out$mods$parms$tc.full,learning.rate=out$mods$lr.mod$lr,bag.fraction=bag.fraction,site.weights=rep(1,times=nrow(dat[c(Split,rep(i,times=sum(dat$response>0)))==i,])),
    #                                   autostop=T,debug.mode=F,silent=!debug.mode,plot.main=F,superfast=F)
    #   #
    #   y <- gbm.interactions(final.mod[[i]])
    #   int <- y$rank.list;
    #   int<-int[int$p<.05,]
    #   int <- int[order(int$p),]
    #   int$p <- round(int$p,4)
    #   names(int) <- c("v1","name1","v2","name2","int.size","p-value")
    #   row.names(int)<-NULL
    #   if(full.fit){
    #     if(nrow(int)>0) out$mods$interactions[[i]] <- int else out$mods$interactions[[i]] <- NULL
    #   }
    # }
    # 
    # if(full.fit) {
    # 
    #   #post processing steps
    #   out$mods$final.mod<-final.mod
    #   var.name<-unlist(lapply(final.mod,function(lst){as.character(lst$contributions[,1])}))
    #   var.contrib<-unlist(lapply(final.mod,function(lst){lst$contributions[,2]}))
    #   var.final<-unique(var.name)
    #   #storing number of variables in final model
    #   out$mods$vnames<- unique(var.name)
    #   #can't take mean here because we need to account for when the variable didn't show up in the model
    #   out$mods$summary<-aggregate(var.contrib,list(Var=var.name),FUN=sum)
    #   out$mods$summary[,2]<-out$mods$summary[,2]/num.splits
    #   names(out$mods$summary)[2]<-"rel.inf"
    #   out$mods$summary<-out$mods$summary[order(out$mods$summary$rel.inf,decreasing=TRUE),]
    #   out$mods$n.vars.final<-length(var.final)
    # 
    #   if(!is.null(unlist(lapply(out$mods$interactions,is.null)))){
    #     interaction.lst<-out$mods$interactions[!unlist(lapply(out$mods$interactions,is.null))]
    #     interactions<-(do.call("rbind",out$mods$interactions))[,1:4] #can't consider p-value just if they were included at least once
    #     out$mods$interactions<-interactions[!duplicated(interactions[,1:4],MARGIN=1),]
    #   } else out$mods$interactions=NULL
    #   return(out)
    # }
    # else return(final.mod)
  #================================================================
  #                        GAM
  #================================================================= 
  
  if(out$modType == "gam") {
    
    # calculating the case weights (equal weights)
    # the order of weights should be the same as presences and backgrounds in the training data
    prNum <- as.numeric(table(dat$Response)["1"])                                # number of presences
    bgNum <- as.numeric(table(dat$Response)["0"])                                # number of backgrounds
    wt <- ifelse(dat$Response == 1, 1, prNum / bgNum)
    
    factor.mask <- na.omit(match(out$factorInputVars,sanitizedVarNames))
    cont.mask <- seq(1:length(sanitizedVarNames))
    if(length(factor.mask)!=0){cont.mask<-cont.mask[-c(factor.mask)]}
 
    if(out$modOptions$AllowShrinkageSmoothers){  
      
      if(out$modOptions$ConsiderLinearTerms){ # creates formula with smooth and linear terms
        startModel = as.formula(paste("Response","~",paste(paste(sanitizedVarNames, collapse=" + "), 
                                                           paste0("s(", sanitizedVarNames, ", bs='ts')",collapse=" + "), sep = " + ")))
        } else { # creates formula with smooth terms oly
          startModel = as.formula(paste("Response","~", paste0("s(", sanitizedVarNames, ", bs='ts')",collapse=" + "), sep = ""))
        }
      
      } else {
        
        if(out$modOptions$ConsiderLinearTerms){ # creates full scope with smooth and linear terms
          startModel = as.formula(paste("Response","~",paste(paste(sanitizedVarNames, collapse=" + "), 
                                                             paste0("s(", sanitizedVarNames, ")",collapse=" + "), sep = " + ")))
        } else { # creates full scope with smooth terms
          startModel = as.formula(paste("Response","~", paste0("s(", sanitizedVarNames, ")",collapse=" + "), sep = ""))
        }
      }
    
    modelGAM <- mgcv::gam(formula = startModel,
                          data = dat,
                          family = binomial(link = "logit"),
                          weights = wt,
                          method = "REML")
    
    return(modelGAM)
  }
  
  # #================================================================
  # #          Habitat Suitability Criterion
  # #================================================================= 
  # if(Model=="udc"){
  #   
  #   out$mods$final.mod[[1]]<-read.udc(out$input$udc.file)
  #   return(out)
  # }
  # 
  # 
  # SplitBackground(out,dat)
  # out$dat$ma$train$Split<-c(Split,rep(0,times=sum(dat$response>0)))
  # #================================================================
  # #                        MARS
  # #================================================================= 
  # if(Model=="mars") {
  #   fit_contribs<-list()
  #   mars.model<-list()
  #   for(i in 1:num.splits){
  #     mars.model[[i]]<-earth(response~.,data=dat[c(Split,rep(i,times=sum(dat$response>0)))==i,],degree=mars.degree,penalty=mars.penalty,glm=list(family=model.family))
  #     
  #     #
  #     if(full.fit) {out$mods$final.mod[[i]]<-mars.model[[i]]
  #     fit_contribs[[i]] <- evimp(mars.model[[i]],trim=TRUE)
  #     }
  #   }
  #   if(full.fit){
  #     assign("fit_contribs",fit_contribs,envir=parent.frame())
  #     return(out)
  #   } else return(mars.model)
  # }

}

### Read Maxent ----------------------------------------------------------------

read.maxent<-function(lambdas){
  lambdas <- read.csv(lambdas,header=FALSE)
  normalizers<-lambdas[(nrow(lambdas)-3):nrow(lambdas),]
  entropy<-normalizers[4,2]
  lambdas<-lambdas[1:(nrow(lambdas)-4),]
  fctType <- rep("raw",times=nrow(lambdas))
  fctType[grep("`",as.character(lambdas[,1]))] <- "reverse.hinge"
  fctType[grep("'",as.character(lambdas[,1]))] <- "forward.hinge"
  fctType[grep("\\^",as.character(lambdas[,1]))]<-"quadratic"
  fctType[grep("[*]",as.character(lambdas[,1]))]<-"product"
  fctType[grep("[(]",as.character(lambdas[,1]))]<-"threshold"
  
  #make these all default to NULL in case the feature type was turned off
  Raw.coef<-Quad.coef<-Prod.coef<-Fwd.Hinge<-Rev.Hinge<-Thresh.val<-Raw.mult<-Quad.mult<-
    Prod.mult<-FH.mult<-FH.cnst<-Rev.mult<-Rev.cnst<-Thresh.cnst<-NULL
  
  if(any(fctType=="raw")){ 
    "Raw.coef"<-lambdas[fctType=="raw",]
    Raw.mult<-c(-sum(Raw.coef[,2]*Raw.coef[,3]/(Raw.coef[,4]-Raw.coef[,3])), Raw.coef[,2]/(Raw.coef[,4]-Raw.coef[,3]))
    Raw.mult[is.nan(Raw.mult)]<-0
  }
  if(any(fctType=="quadratic")){
    "Quad.coef"<-lambdas[fctType=="quadratic",]
    Quad.mult<-c(-sum(Quad.coef[,2]*Quad.coef[,3]/(Quad.coef[,4]-Quad.coef[,3])), Quad.coef[,2]/(Quad.coef[,4]-Quad.coef[,3]))
  }
  if(any(fctType=="product")){ 
    "Prod.coef"<-lambdas[fctType=="product",]
    Prod.coef[,1]<-gsub("[*]",":",Prod.coef[,1])
    Prod.mult<-c(-sum(Prod.coef[,2]*Prod.coef[,3]/(Prod.coef[,4]-Prod.coef[,3])), Prod.coef[,2]/(Prod.coef[,4]-Prod.coef[,3]))
  }
  if(any(fctType=="forward.hinge")){ 
    "Fwd.Hinge"<-lambdas[fctType=="forward.hinge",]
    Fwd.Hinge[,1]<-gsub("'","",Fwd.Hinge[,1])
    FH.mult<-Fwd.Hinge[,2]/(Fwd.Hinge[,4]-Fwd.Hinge[,3])
    FH.cnst<- -Fwd.Hinge[,2]*Fwd.Hinge[,3]/(Fwd.Hinge[,4]-Fwd.Hinge[,3])
  }
  if(any(fctType=="reverse.hinge")){ 
    "Rev.Hinge"<-lambdas[fctType=="reverse.hinge",]
    Rev.Hinge[,1]<-gsub("`","",Rev.Hinge[,1])
    Rev.mult<- -Rev.Hinge[,2]/(Rev.Hinge[,4]-Rev.Hinge[,3])
    Rev.cnst<-Rev.Hinge[,2]*Rev.Hinge[,4]/(Rev.Hinge[,4]-Rev.Hinge[,3])
  }
  if(any(fctType=="threshold")){ 
    "Thresh.val"<-lambdas[fctType=="threshold",]
    Thresh.cnst<-Thresh.val[,2]
  }
  
  
  retn.lst<-list(Raw.coef=Raw.coef,Quad.coef=Quad.coef,Prod.coef=Prod.coef,Fwd.Hinge=Fwd.Hinge,Rev.Hinge=Rev.Hinge,Thresh.val=Thresh.val,Raw.mult=Raw.mult,Quad.mult=Quad.mult,
                 Prod.mult=Prod.mult,FH.mult=FH.mult,FH.cnst=FH.cnst,Rev.mult=Rev.mult,Rev.cnst=Rev.cnst,Thresh.cnst=Thresh.cnst,normalizers=normalizers,entropy=entropy)
  return(retn.lst)
}

### Estimate Learning Rate [for BRT] -------------------------------------------
  # this function estimates optimal number of trees at a variety of learning rates 
  # the learning rate that produces closest to 1000 trees is selected 

est.lr <- function(dat, out){
  
  # set n-folds
  if(out$validationOptions$CrossValidate){
    nFolds <- out$validationOptions$NumberOfFolds
  } else { nFolds <- 10 }
  
  # calculating the case weights
  prNum <- as.numeric(table(dat$Response)["1"])                                 # number of presences
  bgNum <- as.numeric(table(dat$Response)["0"])                                 # number of backgrounds
  wt <- ifelse(dat$Response == 1, 1, prNum / bgNum)

  n.trees <- c(50,100,200,400,800,900,1000,1100,1200,1500,1800,2400,3200)
  lrs <- c(.05,.02,.01,.005,.0025,.001,.0005,.0001) # .1
  lr.out <- NULL
  trees.fit <- 0
  i <- 1
  
  while(trees.fit < 1000 & i <= length(lrs)){
    n <- 1
    repeat {
      try(
        gbm.fit <- gbm.step(data = dat,
                          gbm.x = 8:ncol(dat),                                  # column indices for covariates
                          gbm.y = 4,                                            # column index for response
                          family = out$modelFamily,
                          tree.complexity = ifelse(prNum < 50, 1, 5),
                          learning.rate = lrs[i],
                          bag.fraction = out$modOptions$BagFraction,
                          max.trees = out$modOptions$MaximumTrees,
                          n.trees = n.trees[n],
                          n.folds = nFolds,                                     # number of cross-validation folds
                          site.weights = wt,
                          plot.main = FALSE)
      )
      if(!is.null(gbm.fit)){
        trees.fit <- gbm.fit$n.trees
        row_i <- cbind(lrs = lrs[i], 
                       n.trees = n.trees[n], 
                       trees.fit = trees.fit,
                       cv.dev = gbm.fit$cv.statistics$deviance.mean)
        lr.out <- rbind(lr.out, row_i)
        if(trees.fit >= 1000) break
      } 
      n <- n+1
      if(n > length(n.trees)) break
    }
    i <- i+1
  }
  
  # pick lr and n.trees that gives closest to 1000 trees in final model 
  if (is.null(lr.out)){
    return(lr.out)
  } else{
    lr.out <- as.data.frame(lr.out)
    ab <- coef(lm(trees.fit~log(lrs), data=lr.out))
    lr <- round(as.numeric(exp((1000-ab[1])/ab[2])),6)
    lr.out$abs <- abs(lr.out$trees.fit-1000)
    lr.out$d.lr <- abs(lr.out$lrs-lr)
    lr.out <- lr.out[order(lr.out$abs,lr.out$d.lr),]
    lr.out <- lr.out[1,]
    return(lr.out)
  } 
}

# MODEL SELECTION AND VALIDATION FUNCTIONS -------------------------------------

## Run Cross Validation --------------------------------------------------------

cv.fct <- function(out,         # out list 
                   nfolds,      # number of cross validation folds
                   # sp.no = 1, 
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
  
  # full (training) data 
  family <- out$modelFamily
  data <- out$data$train
  preds <- out$data$train$predicted
  if(out$modType == "brt"){data$predicted <- NULL}
  xdat <- subset(data, select = out$inputVars)
  obs <- out$data$train$Response
  site.weights <- out$data$train$Weight
  selector <- data$ModelSelectionSplit
  n.cases <- nrow(data)
  
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
  nk <- nfolds 
  
  fitted.values <- rep(0, n.cases)
  cor.mat <- matrix(nrow = length(out$inputVars), ncol=nk)
  rownames(cor.mat) <- out$inputVars
  thresh <- NULL
  
  subset.test <- rep(0,nk)
  subset.calib <- as.data.frame(matrix(0,ncol=5,nrow=nk))
  names(subset.calib) <- c("intercept","slope","test1","test2","test3")
  subset.resid.deviance <- rep(0,nk)
  
  ## Start cross validation Fold loop
  
  for (i in 1:nk) {
   
    model.mask <- selector != i  #used to fit model on majority of data
    pred.mask <- selector == i   #used to identify the with-held subset
    # assign("species.subset", obs[model.mask], pos = 1)
    # assign("predictor.subset", xdat[model.mask, ], pos = 1)
    
    if(out$modType == "maxent"){ # prepare input files
      
      if(!file.exists(file.path(out$tempDir, "CVsplits"))){dir.create(file.path(out$tempDir, "CVsplits"))}
      
      data[model.mask,] %>%
        mutate(Species = case_when(Response == 1 ~ "species")) %>%
        drop_na(Species) %>%
        select(-SiteID, -Response, -UseInModelEvaluation, -ModelSelectionSplit, -Weight, -predicted) %>%
        relocate(Species, .before = X) %>%
        write.csv(file.path(out$tempDir, "CVsplits", "training-swd.csv"), row.names = F)
      
      # out$data$background %>%
      #   filter(ModelSelectionSplit != i) %>%
      data[model.mask,] %>%
        mutate(Species = case_when(Response != 1 ~ "background")) %>%
        drop_na(Species) %>%
        select(-SiteID, -Response, -UseInModelEvaluation, -ModelSelectionSplit, -Weight, -predicted) %>%
        relocate(Species, .before = X) %>%
        write.csv(file.path(out$tempDir, "CVsplits", "background-swd.csv"), row.names = F)
    }
    
    # fit new model 
    cv.final.mod <- fitModel(dat = data[model.mask,],
                             out = out,
                             weight = site.weights[model.mask],
                             fullFit = F) 
    
    
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
                                               opt.methods = out$modOptions$thresholdOptimization)[2])
    
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
  
  data$predicted <- fitted.values
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
                  cor.mat = cor.mat, thresh=thresh, cv.train.data = data) # , resp.curves=resp.curves,
  out$cvResults <- cv.list
  return(out)
}


### Calculate Deviance function -----------[see helper functions]---------------

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
    stop("ROC function failed - observed and predicted response must be equal lengths")
  n.x <- as.numeric(length(obs[obs == 0]))
  n.y <- as.numeric(length(obs[obs == 1]))
  xy <- c(preds[obs == 0], preds[obs == 1])
  rnk <- rank(xy)
  wilc <- ((n.x * n.y) + ((n.x * (n.x + 1))/2) - sum(rnk[1:n.x]))/(n.x * n.y)
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
    updateRunLog(paste0("Range of response variable is ", round(pred.range, 2), ". Check family specification."))
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

## Make Model Evaluation Plots -------------------------------------------------

makeModelEvalPlots <- function(out = out){ # previous function name: make.auc.plot.jpg
  
  standResidualFile <- file.path(out$tempDir, paste0(out$modType, "_StandardResidualPlots.png"))
  variableImportanceFile <- file.path(out$tempDir, paste0(out$modType, "_VariableImportance.png"))
  confusionMatrixFile <- file.path(out$tempDir, paste0(out$modType, "_ConfusionMatrix.png"))
  residualSmoothFile <- file.path(out$tempDir, paste0(out$modType, "_ResidualSmoothPlot.png"))
  residualSmoothFctFile <- file.path(out$tempDir, paste0(out$modType, "_ResidualSmoothFunction.rds"))
  ROCAUCFile <- file.path(out$tempDir, paste0(out$modType, "_ROCAUCPlot.png"))  
  AUCPRFile <- file.path(out$tempDir, paste0(out$modType, "_AUCPRPlot.png"))
  calibrationFile <- file.path(out$tempDir, paste0(out$modType, "_CalibrationPlot.png"))
  residualPlotFile <- file.path(out$tempDir, paste0(out$modType, "_PoissonResidualPlots.png"))

  
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
                                                     opt.methods = out$modOptions$thresholdOptimization)[2])
    if(!is.null(out$data$test)){
      
      out$testThresh <- as.numeric(optimal.thresholds(data.frame(ID = 1:length(out$data$test$Response),
                                                                 pres.abs = out$data$test$Response,
                                                                 pred = out$data$test$predicted),
                                                      opt.methods = out$modOptions$thresholdOptimization)[2])
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
  Stats <- list()
  
  if(out$validationOptions$CrossValidate){
    for (i in 1:out$validationOptions$NumberOfFolds){
      Stats[[i]] <- calcStat(x = out$cvResults$cv.train.data[out$cvResults$cv.train.data$ModelSelectionSplit == i,],  
                             family = out$modelFamily,
                             thresh = out$cvResults$thresh[i], 
                             has.split = hasSplit) 
    }
    names(Stats) <- c(1:out$validationOptions$NumberOfFolds)
  }
  
  Stats$train <- calcStat(x = out$data$train, family = out$modelFamily, thresh = out$trainThresh, has.split = hasSplit)
  
  if(out$validationOptions$SplitData){
    Stats$test <- calcStat(x = out$data$test, family = out$modelFamily, thresh = out$testThresh, has.split = hasSplit)
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
  
  ### Create residual surface of input data  

    if(splitType %in% c("none", "cv")){ # use training data
      
      dev.contrib <- calc.deviance(obs = out$data$train$Response, 
                                   preds = out$data$train$predicted,
                                   weights = out$data$train$Weight,
                                   family = out$modelFamily,
                                   return.list = T)$dev.cont

      residualSmoothFct <- resid.image(dev.contrib = dev.contrib, 
                                         dat = out$data$train,
                                         file.name = residualSmoothFile,
                                         label = splitType)
      
    } else { # use testing data
      
      dev.contrib <- calc.deviance(obs = out$data$test$Response, 
                                   preds = out$data$test$predicted,
                                   weights = out$data$test$Weight,
                                   family = out$modelFamily,
                                   return.list = T)$dev.cont
      
      residualSmoothFct <- resid.image(dev.contrib = dev.contrib, 
                                       dat = out$data$test,
                                       file.name = residualSmoothFile,
                                       label = "eval")
    }
  
  saveRDS(object = residualSmoothFct, file = residualSmoothFctFile)
  
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
    
    if(splitType %in% c("cv", "cv/train/test")){
      
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
      
    } # Should above code only be created if CV = TRUE ?? (if yes, keep above if statement, otherwise remove brackets)
    
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
  
  if(splitType == "cv/train/test"){
    
    capture.output(cat("\n\n============================================================",
                       "\n\nEvaluation Statistics"), file = file.path(out$tempDir, paste0(out$modType, "_output.txt")), append = TRUE)
    
    capture.stats(Stats[1:(out$validationOptions$NumberOfFolds)], file.name = file.path(out$tempDir, paste0(out$modType, "_output.txt")), 
                  label = "cv", family = out$modelFamily, opt.methods = out$modOptions$thresholdOptimization, out)
    
    capture.output(cat("\n\n============================================================",
                       "\n\nEvaluation Statistics"), file = file.path(out$tempDir, paste0(out$modType, "_output.txt")), append = TRUE)
    
    capture.stats(Stats["test"], file.name = file.path(out$tempDir, paste0(out$modType, "_output.txt")), 
                  label = "train/test", family = out$modelFamily, opt.methods = out$modOptions$thresholdOptimization, out)
  }
  if(splitType == "cv"){
    
    capture.output(cat("\n\n============================================================",
                       "\n\nEvaluation Statistics"), file = file.path(out$tempDir, paste0(out$modType, "_output.txt")), append = TRUE)
    
    capture.stats(Stats[1:(out$validationOptions$NumberOfFolds)], file.name = file.path(out$tempDir, paste0(out$modType, "_output.txt")), 
                  label = splitType, family = out$modelFamily, opt.methods = out$modOptions$thresholdOptimization, out)
  }
  if(splitType == "train/test"){
    
    capture.output(cat("\n\n============================================================",
                       "\n\nEvaluation Statistics"), file = file.path(out$tempDir, paste0(out$modType, "_output.txt")), append = TRUE)
    
    capture.stats(Stats["test"], file.name = file.path(out$tempDir, paste0(out$modType, "_output.txt")), 
                  label = splitType, family = out$modelFamily, opt.methods = out$modOptions$thresholdOptimization, out)
  }  
  if(splitType == "none"){
    
    capture.output(cat("\n\n============================================================",
                       "\n\nEvaluation Statistics"), file = file.path(out$tempDir, paste0(out$modType, "_output.txt")), append = TRUE)
    
    capture.stats(Stats, file.name = file.path(out$tempDir, paste0(out$modType, "_output.txt")), 
                  label = splitType, family = out$modelFamily, opt.methods = out$modOptions$thresholdOptimization, out)
  }
  
  return(out)
}

### Calculate statistics function ----------------------------------------------

calcStat <- function(x,       # x <- out$data[[i]]
                     family,
                     thresh,  # out$trainThresh  or out$cvResults$thresh
                     has.split){
  
  auc.data <- data.frame(ID=1:nrow(x),pres.abs=x$Response, pred=x$predicted)
  
  p.bar <- sum(auc.data$pres.abs * x$Weight) / sum(x$Weight)
  n.pres <- sum(auc.data$pres.abs>=1)
  n.abs <- nrow(auc.data)-n.pres
  if(has.split & !all(x$UseInModelEvaluation)){ # we're in the train split here for pseudoAbs
    dev.vect <-vector()
    null.dev <-vector()
    for(i in 1:(length(unique(x$ModelSelectionSplit))-1)){
      dev.vect[i] <- calc.deviance(c(rep(0,times=sum(x$ModelSelectionSplit==i)),x$Response[x$Response==1]),
                                 c(x$predicted[x$ModelSelectionSplit==i],x$predicted[x$Response==1]))
      null.dev[i] <- calc.deviance(c(rep(0,times=sum(x$ModelSelectionSplit==i)),x$Response[x$Response==1]),
                                 rep(mean(c(rep(1,times=sum(x$Response==1)),rep(0,times=sum(x$UseInModelEvaluation==i)))),
                                     times=(sum(x$Response==1)+sum(x$ModelSelectionSplit==i))))
    }
    null.dev <- mean(null.dev)
    dev.fit <- mean(dev.vect)
  } else if(has.split & all(x$UseInModelEvaluation)){ # in the test split for pseudoAbs -- deviance shouldn't be calculated because calibration is wrong
    null.dev = NA
    dev.fit = NA
  } else{
    null.dev <- calc.deviance(auc.data$pres.abs, rep(p.bar,times=length(auc.data$pres.abs)), x$Weight, family=family) #*nrow(x$dat)
    if(is.nan(null.dev)){ null.dev <- NA }
    
    dev.fit <- calc.deviance(auc.data$pres.abs, auc.data$pred, x$Weight, family=family) #*nrow(x$dat) Elith does not include this it might cause a weighting issue when averaging I'm not sure if I should include it
    if(is.nan(dev.fit)){ dev.fit <- NA }
    }
  
  dev.exp <- null.dev - dev.fit
  pct.dev.exp <- dev.exp/null.dev*100
  correlation <- cor(auc.data$pres.abs, auc.data$pred)
  
  # have to use roc here because auc in the PresenceAbsence package incorrectly assumes that the value must be greater than .5
  # this isn't necessarily true for an independent evaluation set
  
  auc.fit <- roc(auc.data$pres.abs, auc.data$pred)
  
  calibration.stats <- calibration(auc.data$pres.abs, auc.data$pred, family = family)
  
  if(family %in% c("binomial","bernoulli")){
    
    auc.pr <- pr.curve(scores.class0 = auc.data$pred, weights.class0 = auc.data$pres.abs)[[3]]
    
    cmx <- cmx(auc.data,threshold=thresh)
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
                calibration.stats=calibration.stats, thresh=thresh))
  }
  
  if(family == "poisson"){
    prediction.error <- sum((auc.data$pres.abs-auc.data$pred)^2)
    return(list(n.pres=n.pres,n.abs=n.abs,
                null.dev=null.dev, dev.fit=dev.fit, dev.exp=dev.exp, 
                pct.dev.exp=pct.dev.exp, 
                correlation=correlation,
                auc.data=auc.data, auc.fit=auc.fit, # auc.pr=auc.pr,
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
  
  if(out$modType == "rf"){ 
    trainPred <- pred.fct(mod = out$finalMod, x = out$data$train, modType = out$modType)
    auc$train <- roc(out$dat$train$Response,trainPred)
  } else { trainPred <- out$data$train$predicted }
  
  # add 1 to the drop in AUC to ensure it's greater than zero then I can normalize
  
  cor.mat[,1] <- unlist(auc$train) - permute.predict(inputVars = out$inputVars,
                                                     dat = out$data$train, 
                                                     preds = trainPred,
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
  variable.importance.csv <- file.path(ssimTempDir, paste0(out$modType,"_VariableImportance.csv"))  
  write.table(cbind(predictor=out$inputVars,cor.mat),
              file = variable.importance.csv,
              row.names=FALSE,col.names=TRUE,quote=FALSE,sep=",")
  
  
  # if cross validation we need to avg across folds otherwise plot for each
  # order by the best in the train split
  
  xright <- as.matrix(cor.mat[order(cor.mat[,"train"],decreasing=FALSE),])
  colnames(xright) <- colnames(cor.mat)
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

### Residual Image function ----------------------------------------------------

resid.image <- function(dev.contrib,
                        dat, 
                        file.name,
                        label,
                        create.image = T){
  
  #produces a map of deviance residuals unless we're using independent evaluation data in which case
  #it produces a map of predicted-observed data with a lowess smooth so the relationship is easier to see
  #if there are more than 2000 points a random sample is drawn to speed up the calculations
  #if weights are available these are used for the lowess surface

  if(length(dat$predicted)>2000){
    samp<-seq(1:length(dat$predicted))[order(runif(length(dat$predicted)))][1:2000]
    dev.contrib<-dev.contrib[samp]
    dat<-dat[samp,]
  }
  
  x <- dat$X
  y <- dat$Y
  wgt <- dat$Weight
  
  # dev.contrib is negative for binomial and bernoulli and positive for poisson
  if(label!="eval"){ z <- sign(dat$Response - dat$predicted)*abs(dev.contrib) 
  } else { z <- dat$Response - dat$predicted }
  
  # Remove all infinite values 
  wgt <- wgt[is.finite(z)]
  x <- x[is.finite(z)]
  y <- y[is.finite(z)]
  z <- z[is.finite(z)]
  
  MinCol <- min(z)
  MaxCol <- max(z)
  col.i <- beachcolours(heightrange=c(min(z),max(z)),sealevel=0,s=1,ncolours=(length(table(z))+1))
  f <- function(a,b){ sqrt((a-b)^2)}
  s1 <- seq(from=MinCol,to=MaxCol,length=length(table(z)))
  col.ind <- apply((outer(s1,z,f)),2,which.min)
  
  a <- loess(z~x*y,weights=wgt, 
             control = loess.control(surface = "direct"))
  x.lim <- rep(seq(from=min(x),to=max(x),length=100),each=100)
  y.lim <- rep(seq(from=min(y),to=max(y),length=100),times=100)
  z <- predict(a,newdata=cbind("x"=x.lim,"y"=y.lim))
  x.lim <- seq(from=min(x),to=max(x),length=100)
  y.lim <- seq(from=min(y),to=max(y),length=100)
  z <- matrix(data=z,ncol=100,nrow=100,byrow=TRUE)
  
  # Plot residual smooth with signed and sized residuals on top
  if(create.image){
    png(file=file.name,width=1000,height=1000,pointsize=13)
    par(oma=c(3,3,3,3))
    layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1), heights=c(1,1))
    image(z,x=x.lim,y=y.lim,
          col=beachcolours(heightrange=c(min(z),max(z)),
                           sealevel=0,s=.5,ncolours=length(table(z))),
          main=paste("Spatial pattern of", ifelse(label!="eval"," deviance residuals\n(magnitude and sign)"," prediction error")),
          xlab="X coordinate",ylab="Y coordinate",cex.main=2.2,cex.axis=1.6,cex.lab=1.8)
    
    points(x,y,bg=col.i[col.ind], pch=21,cex=abs(dev.contrib)*1.5)
    par(mar = c(3,2.5,2.5,2))
    
    colrange<-seq(from=MinCol,to=MaxCol,length=100)
    image(1,colrange,
          matrix(data=colrange, ncol=length(colrange),nrow=1),
          col=beachcolours(heightrange=c(MinCol,MaxCol),sealevel=0,ncolours=length(colrange)),
          xlab="",ylab="",
          xaxt="n",cex.main=2,cex.axis=2,cex.lab=2)
    graphics.off()
  }
  return(a)
}

#### beach colours function -----------------------------------------------------

beachcolours <- function (heightrange, 
                          sealevel = 0, 
                          monochrome = FALSE, 
                          s=1,
                          ncolours = if (monochrome) 16 else 64){
  
  if (monochrome)
    return(grey(seq(0, 1, length = ncolours)))
  stopifnot(is.numeric(heightrange) && length(heightrange) == 2)
  stopifnot(all(is.finite(heightrange)))
  depths <- heightrange[1]
  peaks <- heightrange[2]
  dv <- diff(heightrange)/(ncolours - 1)
  epsilon <- dv/2
  lowtide <- max(sealevel - epsilon, depths)
  hightide <- min(sealevel + epsilon, peaks)
  countbetween <- function(a, b, delta) {
    max(0, round((b - a)/delta))
  }
  nsea <- countbetween(depths, lowtide, dv)
  nbeach <- countbetween(lowtide, hightide, dv)
  nland <- countbetween(hightide, peaks, dv)
  colours <- character(0)
  if (nsea > 0)
    colours <- rev(rainbow(nsea, start = 3/6, end = 4/6,s=s))
  if (nbeach > 0)
    colours <- c(colours, rev(rainbow(nbeach, start = 3/12,
                                      end = 5/12,s=s)))
  if (nland > 0)
    colours <- c(colours, rev(rainbow(nland, start = 0, end = 1/6,s=s)))
  return(colours)
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
    } else {stop("Test/Train ROC failed - 'dat' must be either data frame or matrix")}}
  
  obs <- dat[, 2]
  if (length(obs[obs == 0]) == 0) { 
    stop("Test/Train ROC failed - no observed absences in dataset, therefore specificity does not",
         "exist, and modeling, much less Area Under the Curve, is not very meaningful")}
  
  if (length(obs[obs == 1]) == 0) {
    stop("Test/Train ROC failed - no observed presences in dataset, therefore sensitivity does not",
         "exist, and modeling, much less Area Under the Curve, is not very",
         "meaningful")}
  
  if (is.logical(find.auc) == FALSE) { 
    stop("Test/Train ROC failed - 'find.auc' must be of logical type")}
  
  if (is.logical(na.rm) == FALSE) {
    stop("Test/Train ROC failed - 'na.rm' must be of logical type")}
  
  if (is.logical(mark.numbers) == FALSE) {
    stop("Test/Train ROC failed - 'mark.numbers' must be of logical type!")}
  
  if (is.logical(add.legend) == FALSE) {
    stop("Test/Train ROC failed - 'add.legend' must be of logical type!")}
  
  if (is.logical(add.opt.legend) == FALSE) {
    stop("Test/Train ROC failed - 'add.opt.legend' must be of logical type")}
  
  if (is.logical(counter.diagonal) == FALSE) {
    stop("Test/Train ROC failed - 'counter.diagonal' must be of logical type")}
  
  if (length(smoothing) != 1) {
    stop("Test/Train ROC failed - 'smoothing' must be a single number greater than or equal to 1")
  } else {
    if (is.numeric(smoothing) == FALSE) {
      stop("Test/Train ROC failed - 'smoothing' must be a single number greater than or equal to 1")
    } else {
      if (smoothing < 1) {
        stop("Test/Train ROC failed - 'smoothing' must be a single number greater than or equal to 1")
      }}}
  
  if (sum(is.na(dat)) > 0) {
    if (na.rm == TRUE) {
      NA.rows <- apply(is.na(dat), 1, sum)
      updateRunLog(paste0("Warning - Test/Train ROC: ", length(NA.rows[NA.rows > 0]), " rows of data were ignored due to NA values"))
      dat <- dat[NA.rows == 0, ]
    } else { return(NA) }}
  
  dat[dat[, 2] > 0, 2] <- 1
  N.models <- ncol(dat) - 2
  if (is.null(obs.prev) == TRUE) {
    obs.prev <- sum(dat[, 2])/nrow(dat)
  }
  
  if (obs.prev < 0 || obs.prev > 1) {
    stop("Test/Train ROC failed - 'obs.prev' must be a number between zero and one")} 
  
  if (obs.prev == 0) {
    updateRunLog("Warning - Test/Train ROC: Because the observed prevalence was zero, results may be strange")}
  
  if (obs.prev == 1) {
    updateRunLog("Warning - Test/Train ROC: Because the observed prevalence was one, results may be strange")}
  
  if (min(which.model) < 1 || sum(round(which.model) != which.model) != 0){
    stop("Test/Train ROC failed - values in 'which.model' must be positive integers")}
  
  if (max(which.model) > N.models) {
    stop("Test/Train ROC failed - values in 'which.model' must not be greater than number of models in 'dat'!")}
  
  if (is.null(model.names) == TRUE) {
    model.names <- if (is.null(names(dat)) == FALSE) {
      names(dat)[-c(1, 2)]
    } else {
      paste("Model", 1:N.models)
    }}
  
  if (N.models != length(model.names) && (length(which.model) != 1 || length(model.names) != 1)) {
    stop("Test/Train ROC failed - If 'model.names' is specified it must either be a single name, or a vector",
         "of the same length as the number of model predictions in 'dat'")}
  
  if (is.null(legend.text) == TRUE) { legend.text <- model.names }
  
  if (length(legend.text) != N.models) {
    stop("Test/Train ROC failed - 'opt.legend.text' must be of same length as 'opt.methods'")}
  
  dat <- dat[, c(1, 2, which.model + 2)]
  if (length(model.names) != 1) { model.names <- model.names[which.model] }
  if (length(legend.text) != 1) { legend.text <- legend.text[which.model] }
  
  N.dat <- ncol(dat) - 2
  if (is.null(obs.prev) == TRUE) { obs.prev <- sum(dat[, 2])/nrow(dat) }
  if (obs.prev < 0 || obs.prev > 1) {
    stop("Test/Train ROC failed - 'obs.prev' must be a number between zero and one")}
  
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
          stop("Test/Train ROC failed - invalid optimization method")
        } else {
          opt.methods <- POSSIBLE.meth[opt.methods]
        }}
      
      if (sum(opt.methods %in% POSSIBLE.meth) != N.meth) {
        stop("Test/Train ROC failed - invalid optimization method") }
      
      if (is.null(opt.legend.text) == TRUE) { opt.legend.text <- opt.methods }
      
      if (length(opt.legend.text) != N.meth) {
        stop("Test/Train ROC failed - 'opt.legend.text' must be of same length as 'opt.methods'") }
      
      if ("ReqSens" %in% opt.methods) {
        if (missing(req.sens)) {
          updateRunLog("Warning - Test/Train ROC: required sensetivity defaults to 0.85")
          req.sens <- 0.85
        }}
      
      if ("ReqSpec" %in% opt.methods) {
        if (missing(req.spec)) {
          updateRunLog("Warning - Test/Train ROC: required specificity defaults to 0.85")
          req.spec <- 0.85
        }}
      
      if ("Cost" %in% opt.methods) {
        if (missing(FPC) || missing(FNC)) {
          updateRunLog("Warning - Test/Train ROC: costs assumed to be equal")
          FPC <- 1
          FNC <- 1
        }
        if (FPC <= 0 || FNC <= 0) { stop("costs must be positive") }
        if (is.logical(cost.line) == FALSE) {
          stop("Test/Train ROC failed - 'cost.line' must be of logical type")
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
      stop("Test/Train ROC failed - 'opt.thresholds' must be 'TRUE', 'FALSE', or numeric") }
    
    if (min(opt.thresholds) < 0) { 
      stop("Test/Train ROC failed - 'opt.thresholds' can not be negative") }
    
    if (max(opt.thresholds) > 1) {
      if (N.thr == 1 && round(opt.thresholds) == opt.thresholds) {
        opt.thresholds <- seq(length = opt.thresholds,
                              from = 0, to = 1)
        N.thr <- length(opt.thresholds)
      } else {
        stop("Test/Train ROC failed - non-interger, non-logical 'opt.thresholds' greater than 1")
      }}
    
    N.opt.thresh <- length(opt.thresholds)
    if (is.null(opt.legend.text)) {
      opt.legend.text <- rep("threshold", N.opt.thresh)
    }
    
    if (length(opt.legend.text) != N.opt.thresh) {
      stop("Test/Train ROC failed - length of 'opt.legend.text' does not match number of specified thresholds") }
    
    if (is.null(pch)) { pch <- 1:N.opt.thresh }
    
    if (length(pch) != N.opt.thresh) {
      stop("Test/Train ROC failed - length of 'pch' does not match number of specified thresholds") }
    
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
  
  # if("try-error" %in% class(gam1) | !gam1$converged){
  #   gam1 <- try(glm(obs ~ ns(preds, df=3), 
  #                   weights=rep(1, length(preds)), 
  #                   family=binomial),silent=TRUE)
  # }
  
  if("try-error" %in% class(gam1)){  # | !gam1$converged
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
  if(label == "none"){ label = "Training" }
  if (label == "cv") {label = "Cross validation"}
  
  capture.output(cat(" applied to",label, "split:\n",sep=" "), file=file.name, append=TRUE)
  
  capture.output(cat( "\n", "\n\t Correlation Coefficient      :",mean(unlist(lapply(Stats.lst,function(lst){lst$correlation}))),
                      if(label == "Cross validation"){
                        paste(" (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$correlation}))),digits=5), ")",sep="")},
                      if(!out$hasSplit |(out$hasSplit & label == "Training")){ 
                        paste("\n\t NULL Deviance                : ", signif(mean(unlist(lapply(Stats.lst,function(lst){lst$null.dev}))),digits=5),
                              if(out$hasSplit & label == "Training"){" (Averaged over background splits)"},
                              if(label == "Cross validation"){
                                paste(" (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$null.dev}))),digits=5), ")",sep="")},
                              "\n\t Fit Deviance                 : ", signif(mean(unlist(lapply(Stats.lst,function(lst){lst$dev.fit}))),digits=5),
                              if(out$hasSplit & label == "Training"){" (Averaged over background splits)"},  
                              if(label == "Cross validation"){
                                paste(" (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$dev.fit}))),digits=5),")",sep="")},
                              "\n\t Explained Deviance           : ", signif(mean(unlist(lapply(Stats.lst,function(lst){lst$dev.exp}))),digits=5),
                              if(label == "Cross validation"){
                                paste(" (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$dev.exp}))),digits=5),")",sep="")},
                              "\n\t Percent Deviance Explained   : ", signif(mean(unlist(lapply(Stats.lst,function(lst){lst$pct.dev.exp}))),digits=5),
                              if(label == "Cross validation"){
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
      if (label %in% c("Training", "Final evaluation")){
        paste("\n\t Threshold                    : ", Stats.lst[[1]]$thresh)
      } else { 
        paste("\n\t Mean Threshold               : ", mean(unlist(lapply(Stats.lst,function(lst){lst$thresh}))),
              " (sd ", signif(sd(unlist(lapply(Stats.lst,function(lst){lst$thresh}))),digits=5),")",sep="")
      },
      "\n\n\t Confusion Matrix: \n\n"),
      if(label %in% c("Training", "Final evaluation")){ print.table(Stats.lst[[1]]$Cmx)
      } else {
        a <-lapply(Stats.lst,function(lst){lst$Cmx})
        cmx <-a[[1]]
        for(i in 2:length(a)){ cmx<-cmx+a[[i]] }
        print.table(cmx)
      },
      cat(
        "\n\t AUC                          : ",mean(unlist(lapply(Stats.lst,function(lst){lst$auc.fit}))),
        if(label == "Cross validation"){
          paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$auc.fit}))),digits=5),")",sep="")},
        "\n\t AUC-pr                       : ",mean(unlist(lapply(Stats.lst,function(lst){lst$auc.pr}))),
        if(label == "Cross validation"){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$auc.pr}))),digits=5),")",sep="")},
        
        "\n\t Percent Correctly Classified : ",mean(unlist(lapply(Stats.lst,function(lst){lst$Pcc}))),
        if(label == "Cross validation"){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$Pcc}))),digits=5),")",sep="")},
        
        "\n\t Sensitivity                  : ",mean(unlist(lapply(Stats.lst,function(lst){lst$Sens}))),
        if(label == "Cross validation"){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$Sens}))),digits=5),")",sep="")},
        
        "\n\t Specificity                  : ",mean(unlist(lapply(Stats.lst,function(lst){lst$Specf}))),
        if(label== "Cross validation"){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$Specf}))),digits=5),")",sep="")},
        
        "\n\t Kappa                        : ",mean(unlist(lapply(Stats.lst,function(lst){lst$Kappa}))),
        if(label == "Cross validation"){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$Kappa}))),digits=5),")",sep="")},
        
        "\n\t True Skill Statistic         : ",mean(unlist(lapply(Stats.lst,function(lst){lst$Tss}))),
        if(label == "Cross validation"){paste(" (sd ",signif(sd(unlist(lapply(Stats.lst,function(lst){lst$Tss}))),digits=5),")",sep="")},
        "\n"),
      
      file=file.path(out$tempDir, paste0(out$modType, "_output.txt")),append=TRUE)
  }
  
  if(!out$pseudoAbs){
    capture.output(cat( "\n\n   Calibration Statistics",
                        "\n\t Intercept (general calibration)                            : ",mean(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[1]}))),
                        if(label == "Cross validation"){paste(" (sd ",
                                                              signif(sd(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[1]}))),digits=5), ")",sep="")},
                        "\n\t Slope   (direction and variation in fit)                   : ",mean(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[2]}))),
                        if(label == "Cross validation"){paste(" (sd ",
                                                              signif(sd(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[2]}))),digits=5),
                                                              ")",sep="")},
                        "\n\t Testa0b1 (overall reliability of predictors)               : ",mean(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[3]}))),
                        if(label == "Cross validation"){paste(" (sd ",
                                                              signif(sd(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[3]}))),digits=5),
                                                              ")",sep="")},
                        "\n\t Testa0|b1(incorrect calibration given correct refinement)  : ",mean(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[4]}))),
                        if(label == "Cross validation"){paste(" (sd ",
                                                              signif(sd(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[4]}))),digits=5),
                                                              ")",sep="")},
                        "\n\t Testb1|a (refinement given correct calibration)            : ",mean(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[5]}))),
                        if(label == "Cross validation"){paste(" (sd ",
                                                              signif(sd(unlist(lapply(Stats.lst,function(lst){lst$calibration.stats[5]}))),digits=5),
                                                              ")",sep="")},
                        
                        "\n\n",
                        file=file.path(out$tempDir, paste0(out$modType, "_output.txt")),append=TRUE))
    #if(label=="crossValidation"){cat("\n\n   Pooled Calibration Statistics\n",print.table(cbind(names(out$cv$pooled.calib),out$cv$pooled.calib)))}
    #something I should include later
  }
}

## response curves function -----------------------------------------------------

response.curves <- function(out){
  
  # Desanitize output variable names
  for (i in 1:length(out$finalVars)) {
    if (grepl("\\s", out$finalVars[i])){
      oldVarName <- out$finalVars[i]
      newVarName <- gsub("`", "", oldVarName)
      out$finalVars[i] <- newVarName
    }
  }
  # if(out$modType %in% c("mars")){ nVars <- nrow(out$mod$summary)
  # if(out$modType %in% c("udc")) nVars <- out$mods$n.vars.final
  if(out$modType %in% c("glm", "rf", "maxent", "brt", "gam")){ nVars <- out$nVarsFinal  }
  
  pcol <- ceiling(sqrt(nVars))
  prow <- ceiling(nVars/pcol)
  
  
  # rf partial plot does something odd with the y axis
  dat <- data.frame(out$data$train[,out$finalVars])
  names(dat) <- out$finalVars
  resp <- out$data$train$Response
  Xp <- as.data.frame(matrix(NA, nc = ncol(dat), nr = nrow(dat),
                             dimnames = list(NULL, colnames(dat))))
  for (i in 1:ncol(dat)) {
    if (is.numeric(dat[, i])) {
      Xp[, i] <- mean(dat[, i])
    }
    else {
      Xp[, i] <- as.factor(rep(names(which.max(summary(dat[,i]))), nrow(dat)))
      levels(Xp[, i]) <- levels(dat[, i])
    }
  }
  
  # dir.create(paste0(out$tempDir, out$modType, "_ResponseCurves"))
  
  # for (k in c(1,2)){
  k <- 1
    if(k==1){ png(file.path(out$tempDir, paste0(out$modType,"_ResponseCurves.png")), width=2000, height=2000, pointsize = 20) 
      par(oma=c(1,6,4,2),mfrow=c(prow,pcol))}                   
    for (i in sort(match(out$finalVars,names(dat)))) {
      if (k==2){ png(filename=file.path(out$tempDir, paste0(out$modType,"_ResponseCurves"), paste0(names(dat)[i],".png",sep="")),width=2000, height=2000, pointsize = 20)
        par(mar=c(5,8,6,1))
      }
      if (!is.factor(dat[, i])) {
        xr <- range(dat[, i])
        Xp1 <- Xp
        Xp1[, i] <- seq(xr[1], xr[2], len = nrow(dat))
        } else {
          Xp1 <- Xp
          Nrepcat <- floor(nrow(dat)/length(levels(dat[,i])))
          Xp1[, i] <- as.factor(c(rep(levels(dat[, i])[1],
                                    nrow(dat) - (Nrepcat * length(levels(dat[,i])))), 
                                rep(levels(dat[, i]), each = Nrepcat)))
      }
      Xf <- matrix(nrow=nrow(Xp1),ncol=1)
      Xf[,1] <- pred.fct(mod = out$finalMod, x = as.data.frame(Xp1), modType = out$modType)
      
      y.lim <- c(0,1)
      if(out$modelFamily == "poisson"){ y.lim <- range(apply(Xf,1,mean)) } 
      
      plot(Xp1[, i],apply(Xf,1,mean), ylim = y.lim, xlab = "",
           ylab = "", type = "l", main = names(dat)[i],lwd=ifelse(k==1,4,6),
           cex.main=3.5,xaxt="n",yaxt="n")
      if(k==2){ mtext("Predicted Value", side=2,line=5,cex=3.5) }
      axis(1,labels=FALSE)
      axis(2,labels=FALSE)
      axis(1,line=0.75,lwd=0,cex.axis=2)
      axis(2,line=0.75,lwd=0,cex.axis=2)
      
      rug(dat[resp==1,i],col="red",lwd=2)
      rug(dat[resp==0,i],col="blue",lwd=2)
      if(k==2){ graphics.off() }  
    } 
    if(k==1){
      mtext("Predicted Value",side=2,cex=3,line=2,outer=TRUE)
      graphics.off()
    }     
  # } # end k statement 
}  
