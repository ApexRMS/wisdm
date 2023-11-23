## ------------------------- 
## wisdm - helper functions  
## ApexRMS, March 2022       
## ------------------------- 


# Calculate Deviance function --------------------------------------------------

calc.deviance <- function(obs,   # observed response
                          preds, # predicted response 
                          weights = rep(1,length(obs)), 
                          family = "binomial", 
                          calc.mean = TRUE, 
                          return.list = FALSE)
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
  
  if (return.list){
    deviance <- list(deviance=deviance,dev.cont=deviance.contribs)
  } 
  return(deviance)
  
}

# Predict Function -------------------------------------------------------------

pred.fct <- function(mod,      # mod = the model fit object
                     x,        # x = data to predict for
                     modType){ # modType = one of mars, glm, rf, brt, maxlike at present 
  
  y <- rep(NA,nrow(x))
  
  if(modType %in% c("glm","mars")){
   y <- try(as.vector(predict(mod, x, type="response")), silent=TRUE) 
  }
   
  if(modType == "rf"){
    # make predictions from complete data only #
    y[complete.cases(x)] <- try(as.vector(predict(mod, newdata=x[complete.cases(x),], type="vote")[,2]), silent=TRUE)
  }
  if(modType == "maxent"){
    y[complete.cases(x)]<-try(maxent.predict(mod, x[complete.cases(x),]), silent=TRUE)
  }
  if(modType == "brt"){
    y[complete.cases(x)] <- try(predict(mod, x[complete.cases(x),], mod$gbm.call$best.trees, type="response"), silent=TRUE)
  }
  if(modType %in% c("gam")){
    y[complete.cases(x)] <- try(predict.gam(mod, x[complete.cases(x)], type="response"), silent=TRUE) 
  }

  # if(modType=="udc"){
  #   
  #   y[complete.cases(x)]<-try(udc.predict(model,x[complete.cases(x),]),silent=TRUE)
  # }
  # if(class(y)=="try-error") stop("Predicting the response for the new values failed.  One probable cause is that you are trying to predict to factor levels that were not present during model fitting.")
  return(y)
} # end pred.vals function 

## glm predict function --------------------------------------------------------

glm.predict <- function(model,x) {
  # retrieve key items from the global environment #
  # make predictions.
  
  y <- as.vector(stats::predict(object = model, newdata = x, type = "response"))
  
  # encode missing values as -1.
  y[is.na(y)]<- NaN
  
  # return predictions.
  return(y)
}

## rf predict function ---------------------------------------------------------

rf.predict <- function(model, x){
  # retrieve key items from the global environment #
  # make predictions from complete data only #
  y <- rep(NA,nrow(x))
  y[complete.cases(x)] <- try(as.vector(predict(mod, newdata=x[complete.cases(x),], type="vote")[,2]), silent=TRUE)
  
  # return predictions.
  return(y)
}

## tweak prediction function (for rf) ------------------------------------------

tweak.p <- function(p){
  p[p==1]<-max(p[p<1])
  p[p==0]<-min(p[p>0])
  return(p)
}

## maxent predict function -----------------------------------------------------

maxent.predict <- function(model,x){
  
  if(names(model[[1]])[1]=="Raw.coef")
  {attach(model[[1]])
    on.exit(detach(model[[1]]))
  } else{attach(model)
    on.exit(detach(model))
  }
  
  # build prediction shell filled with NA
  prediction <- rep(NA, nrow(x))
  names(prediction) <- row.names(x)
  
  # These all default to zero in case the feature was excluded
  Raw <- Quad <- Prod <- Forw <- Rev <- Thresh <- 0
  
  if(!is.null(Raw.coef)){
    Raw<- model.matrix(as.formula(paste("~",paste(Raw.coef[,1],collapse=" + "),sep="")),x)
    Raw<-apply(t(Raw)*Raw.mult,2,sum)
  }
  if(!is.null(Quad.coef)){  
    Quad<- model.matrix(as.formula(paste("~ I(",paste(Quad.coef[,1],collapse=") + I("),")",sep="")),x)
    Quad<-apply(t(Quad)*Quad.mult,2,sum)
  }
  if(!is.null(Prod.coef)){  
    Prod<- model.matrix(as.formula(paste("~ ",paste(Prod.coef[,1],collapse=" + "),sep="")),x)
    Prod<-apply(t(Prod)*Prod.mult,2,sum)
  }
  if(!is.null(Fwd.Hinge)){
    Forw<- model.matrix(as.formula(paste("~ ",paste("-1 + ",paste("I(",paste(Fwd.Hinge[,1],paste(Fwd.Hinge[,1],
                                                                                                 paste("(",Fwd.Hinge[,3],")",sep=""),sep=">"),sep="*("),collapse=")) + "),"))",sep=""),sep="")),x)
    Forw.Cst<-Forw!=0
    Forw<-apply(t(Forw)*FH.mult,2,sum)+apply(t(Forw.Cst)*FH.cnst,2,sum)
  }
  if(!is.null(Rev.Hinge)){
    Rev<- model.matrix(as.formula(paste("~ ",paste("-1 + ",paste("I(",paste(Rev.Hinge[,1],paste(Rev.Hinge[,1],
                                                                                                paste("(",Rev.Hinge[,4],")",sep=""),sep="<"),sep="*("),collapse=")) + "),"))",sep=""),sep="")),x)
    Rev.Cst<-Rev!=0
    Rev <- apply(t(Rev)*Rev.mult,2,sum) + apply(t(Rev.Cst)*Rev.cnst,2,sum)
  }
  if(!is.null(Thresh.val)){
    Thresh.val[,1]<-gsub("=","==",Thresh.val[,1])
    Thresh<-model.matrix(as.formula(paste("~ ",paste("I",Thresh.val[,1],collapse=" + ",sep=""),sep="")),x)
    Thresh<-apply(t(Thresh[,2:ncol(Thresh)])*Thresh.cnst,2,sum)
  }
  
  S <- Raw + Quad + Prod + Forw + Rev + Thresh - normalizers[1,2]
  
  qx <-  exp(S)/normalizers[2,2]
  predVals <- qx*exp(entropy)/(1+qx*exp(entropy))
  prediction[names(predVals)] <- predVals
  return(prediction)
}

## brt predict function --------------------------------------------------------

brt.predict <- function(model,x) {
  
  y <- rep(NA,nrow(x))
  y[complete.cases(x)] <- gbm::predict.gbm(model, x[complete.cases(x),], model$gbm.call$best.trees, type="response")

  # # make predictions from full data #
  # y <- predict.gbm(model, x, model$target.trees, type="response")
  # # encode missing values as -1.
  # a <- complete.cases(x)
  # y[!(a)]<- NaN

  # return predictions.
  return(y)
}

## gam predict function --------------------------------------------------------

gam.predict <- function(model,x) {
  
  y <- rep(NA,nrow(x))
  y[complete.cases(x)] <- mgcv::predict.gam(model, x[complete.cases(x),], type="response")
  
  # return predictions.
  return(y)
}
# ## mars predict function -------------------------------------------------------
# 
# mars.predict <- function(model,x) {
#   # retrieve key items from the global environment #
#   # make predictions.
#   y <- rep(NA,nrow(x))
#   y[complete.cases(x)] <- as.vector(predict.mars(model,x[complete.cases(x),])$prediction[,1])
#   
#   #which(is.na(y)
#   # encode missing values as -1.
#   y[is.na(y)]<- NaN
#   
#   # return predictions.
#   return(y)
# }
# 
# ### predict mars [nested] function ---------------------------------------------
#   
#   predict.mars <- function (mars.glm.object,new.data){
#     
#     # version 3.1 - developed in R 2.3.1 using mda 0.3-1
#     #
#     # calculates a mars/glm object in which basis functions are calculated
#     # using an initial mars model with single or multiple responses
#     # data for individual species are then fitted as glms using the
#     # common set of mars basis functions with results returned as a list
#     #
#     # takes as input a dataset and args selecting x and y variables, and degree of interaction
#     # along with site and species weights, the CV penalty, and the glm family argument
#     # the latter would normally be one of "binomial" or "poisson" - "gaussian" could be used
#     # but in this case the model shouldn't differ from that fitted using mars on its own
#     #
#     # naming problem for dataframes fixed - je - 15/12/06
#     #
#     # requires mda and leathwick/elith's mars.export
#     
#     require(mda)
#     
#     # first recreate both the original mars model and the glm model
#     
#     # setup input data and create original temporary data
#     
#     dataframe.name <- mars.glm.object$mars.call$dataframe  # get the dataframe name
#     mars.x <- mars.glm.object$mars.call$mars.x
#     mars.y <- mars.glm.object$mars.call$mars.y
#     n.spp <- length(mars.y)
#     family <- mars.glm.object$mars.call$family
#     mars.degree <- mars.glm.object$mars.call$degree
#     penalty <- mars.glm.object$mars.call$penalty
#     site.weights <- mars.glm.object$weights[[1]]
#     spp.weights <- mars.glm.object$weights[[2]]
#     
#     print("creating original data frame...",quote=FALSE)
#     
#     base.data <- as.data.frame(eval.parent(parse(text = dataframe.name),n=3)) #aks
#     
#     x.temp <- eval(base.data[, mars.x])                 #form the temporary datasets
#     base.names <- names(x.temp)
#     
#     xdat <- mars.new.dataframe(x.temp)[[1]]
#     
#     ydat <- as.data.frame(base.data[, mars.y])
#     names(ydat) <- names(base.data)[mars.y]
#     
#     assign("xdat", xdat, pos = 1)               #and assign them for later use
#     assign("ydat", ydat, pos = 1)
#     
#     # now create the temporary dataframe for the new data
#     
#     print("checking variable matching with new data",quote = FALSE)
#     
#     new.names <- names(new.data)
#     
#     for (i in 1:length(base.names)) {
#       
#       name <- base.names[i]
#       
#       if (!(name %in% new.names)) {
#         print(paste("Variable ",name," missing from new data",sep=""),quote = FALSE)  #aks
#         return()
#       }
#     }
#     
#     print("and creating temporary dataframe for new data...",quote=FALSE)
#     
#     selector <- match(names(x.temp),names(new.data))
#     
#     pred.dat <- mars.new.dataframe(new.data[,selector])[[1]]
#     
#     assign("pred.dat", pred.dat, pos = 1)               #and assign them for later use
#     
#     # fit the mars model and extract the basis functions
#     
#     print(paste("re-fitting initial mars model for",n.spp,"responses"),quote = FALSE)
#     print(paste("using glm family of",family),quote = FALSE)
#     
#     #mars.fit <- mars(x = xdat, y = ydat, degree = mars.degree, w = site.weights,
#     #  wp = spp.weights, penalty = penalty)
#     
#     mars.fit <- mars.glm.object$mars.object  #AKS
#     
#     old.bf.data <- as.data.frame(eval(mars.fit$x))
#     n.bfs <- ncol(old.bf.data)
#     bf.names <- paste("bf", 1:n.bfs, sep = "")
#     old.bf.data <- as.data.frame(old.bf.data[,-1])
#     names(old.bf.data) <- bf.names[-1]
#     
#     new.bf.data <- as.data.frame(mda:::model.matrix.mars(mars.fit,pred.dat))
#     new.bf.data <- as.data.frame(new.bf.data[,-1])
#     names(new.bf.data) <- bf.names[-1]
#     
#     # now cycle through the species fitting glm models
#     
#     print("fitting glms for individual responses", quote = F)
#     
#     prediction <- as.data.frame(matrix(0, ncol = n.spp, nrow = nrow(pred.dat)))
#     names(prediction) <- names(ydat)
#     standard.errors <- as.data.frame(matrix(0, ncol = n.spp, nrow = nrow(pred.dat)))
#     names(standard.errors) <- names(ydat)
#     
#     for (i in 1:n.spp) {
#       
#       print(names(ydat)[i], quote = FALSE)
#       model.glm <- glm(ydat[, i] ~ ., data = old.bf.data, weights = site.weights,
#                        family = family, maxit = 100)
#       temp <- predict.glm(model.glm,new.bf.data,type="response",se.fit=TRUE)
#       prediction[,i] <- temp[[1]]
#       standard.errors[,i] <- temp[[2]]
#       
#     }
#     
#     return(list("prediction"=prediction,"ses"=standard.errors))
#   }
# 
# #### mars new dataframe [nested] function --------------------------------------
# 
# mars.new.dataframe <- function(input.data){
#     
#     # j leathwick, j elith - August 2006
#     #
#     # version 3.1 - developed in R 2.3.1 using mda 0.3-1
#     #
#     # takes an input data frame and checks for factor variables
#     # converting these to dummy variables, one each for each factor level
#     # returning it for use with mars.glm so that factor vars can be included
#     # in a mars analysis
#     
#     
#     if (!is.data.frame(input.data)) {
#       print("input data must be a dataframe..",quote = FALSE)
#       return()
#     }
#     
#     n <- 1
#     for (i in 1:ncol(input.data)) {  #first transfer the vector variables
#       if (is.vector(input.data[,i])) {
#         if (n == 1) {
#           output.data <- as.data.frame(input.data[,i])
#           names.list <- names(input.data)[i]
#           var.type <- "vector"
#           factor.level <- "na"
#         }
#         else {
#           output.data[,n] <- input.data[,i]
#           names.list <- c(names.list,names(input.data)[i])
#           var.type <- c(var.type,"vector")
#           factor.level <- c(factor.level,"na")
#         }
#         names(output.data)[n] <- names(input.data)[i]
#         n <- n + 1
#       }
#     }
#     
#     for (i in 1:ncol(input.data)) {  # and then the factor variables
#       if (is.factor(input.data[,i])) {
#         temp.table <- summary(input.data[,i])
#         for (j in 1:length(temp.table)) {
#           names.list <- c(names.list,names(input.data)[i])
#           var.type <- c(var.type,"factor")
#           factor.level <- c(factor.level,names(temp.table)[j])
#           output.data[,n] <- ifelse(input.data[,i] == names(temp.table)[j],1,0)
#           names(output.data)[n] <- paste(names(input.data)[i],".",names(temp.table)[j],sep="")
#           n <- n + 1
#         }
#       }
#     }
#     
#     lineage <- data.frame(names(output.data),names.list,var.type,factor.level)
#     for (i in 1:4) lineage[,i] <- as.character(lineage[,i])
#     names(lineage) <- c("full.name","base.name","type","level")
#     
#     return(list(dataframe = output.data, lineage = lineage))
#   }
