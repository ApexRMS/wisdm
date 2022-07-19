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

## glm predict function --------------------------------------------------------

glm.predict <- function(model,x) {
  # retrieve key items from the global environment #
  # make predictions.
  
  y <- as.vector(stats::predict(object = model, newdata = x ,type = "response"))
  
  # encode missing values as -1.
  y[is.na(y)]<- NaN
  
  # return predictions.
  return(y)
}

# ## brt precict function --------------------------------------------------------
# 
# brt.predict <- function(model,x) {
#   # retrieve key items from the global environment #
#   # make predictions from complete data only #
#   #y <- rep(NA,nrow(x))
#   #y[complete.cases(x)] <- predict.gbm(model, x[complete.cases(x),],model$target.trees,type="response")
#   
#   # make predictions from full data #
#   y <- predict.gbm(model,x,model$target.trees,type="response")
#   # encode missing values as -1.
#   a<-complete.cases(x)
#   y[!(a)]<- NaN
#   
#   # return predictions.
#   return(y)
# }
# 
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
