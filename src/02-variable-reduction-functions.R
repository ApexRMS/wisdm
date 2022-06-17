## sdsim - variable reduction functions
## ApexRMS, May 2022

# pairs explore function -------------------------------------------------------

pairsExplore <- function(inputData,     # dataframe with response and covariate values per site 
                         selectedCovs,  # covariates to consider
                         options,       # covariate selection options
                         factorVars,
                         family,
                         outputFile,
                         # pres=TRUE,
                         # absn=TRUE,
                         # bgd=TRUE,
                         # Debug=FALSE,
                         seed=1){
  # input data
  response <- inputData$Response
  dat <- select(inputData, all_of(selectedCovs))
  
  # input options
  displayHighCors <- options$DisplayHighestCorrelations
  minCor <- options$CorrelationThreshold
  numPlots <- options$NumberOfPlots
  
  # calculate the proportion of missing data per variable  
  missing.summary <- 1-apply(apply(dat,2,complete.cases),2,sum)/nrow(dat)
  
  # subsample the data so we can calculate correlations quickly
  response.table <- table(response)
  max.points <- 1500
  if(any(response.table > max.points)){
    for(i in names(response.table[response.table > max.points])){
      s <- sample(which(response==i,arr.ind=TRUE),size=(sum(response==i) - max.points))
      dat <- dat[-c(s),]
      response <- response[-c(s)]
      # TrueResponse <- TrueResponse[-c(s)]
    }
  }   
  #the deviance calculation requires even columns which will be removed for the pairs explore
  #but to get the same answer for the plot I need the same subsample
  
  for.dev <- list(dat=dat,response=response)  
  
  ## Marking factor columns as such after calculating deviance these will be removed
  if(length(factorVars)>0){
    for(i in factorVars) dat[,i] <- factor(dat[,i])
  }
  
  devExp <- vector()
  for(i in (1:ncol(for.dev$dat))){
    devExp[i] <- try(my.panel.smooth(x = for.dev$dat[,i], 
                                     y = for.dev$response,
                                     plot.it=FALSE,
                                     family=family),silent=TRUE)
  }
  devInfo <- as.data.frame(devExp,row.names=names(for.dev[[1]]))
  # write.csv(devInfo, file = paste(dirname(output.file),"devInfo.csv",sep="/"))
  
  # now removing factor columns
  if(length(factorVars)>0){ 
    dat <- select(dat, -factorVars)
    warning("The covariate correlation tool does not consider categorical predictors")
  }
  
  #after calculating the deviance for all predictors we have to remove the excluded predictors for the following plots
  for.dev$dat <- dat 
  
  # Remove columns with only one unique value
  if(is.null(dim(dat))){
    stop("You must have more than 1 non-categorical predictor to use the pairs plot") 
  }
    
  dat <- try(dat[,as.vector(apply(dat,2,var,na.rm=TRUE)==0)!=1],silent=TRUE)
  if(class(dat)=="try-error"){ stop("Site data contains nonnumeric values please remove and continue") }
  
  # record correlations for later plots
  cmat <- cor(dat, use="pairwise.complete.obs")
  smat <- cor(dat, method="spearman", use="pairwise.complete.obs")
  if(dim(dat)[1]<2000){
    kmat <- cor(dat, method="kendall", use="pairwise.complete.obs")
    } else { 
      s <-sample(seq(1:dim(dat)[1]),size=2000,replace=FALSE)
      kmat <- cor(dat[s,], method="kendall", use="pairwise.complete.obs")
      }
  cmat <- pmax(abs(cmat), abs(smat), abs(kmat), na.rm=TRUE)
  cmat[is.na(cmat)] <- 0
  High.cor <- sort(apply(abs(cmat) > minCor,2,sum)-1, decreasing=TRUE)
  
  # take the top num.plots to put in the pairs plot or if the looking at a single
  # predictor and other predictors it's correlated with, take the top num.plots-1
  # of those with which it is correlated
  if(displayHighCors){
    # take the column of the correlation matrix corresponding to the
    # predictor with the highest number of total correlations record the names
    # of the predictors that are correlated with this one predictor
    temp <- cmat[rownames(cmat)==names(High.cor[1]),]
    corWHigh <- temp[abs(cmat[,colnames(cmat) == names(High.cor[1])]) > minCor]
    
    #record counts of total number of correlations with all predictors for those
    #predictors that are highly correlated with the Highest predictor
    # High.cor <- sort(High.cor[names(corWHigh)], decreasing=TRUE)    
  }
  HighToPlot <- dat[,match(names(High.cor),names(dat))[1:min(numPlots,length(High.cor))]]
  for.dev$dat <- for.dev$dat[,match(names(High.cor),names(dat))[1:min(numPlots,length(High.cor))]]
  
  cor.hightoplot <- abs(cor(HighToPlot,use="pairwise.complete.obs"))
  diag(cor.hightoplot) <- 0
  cor.hightoplot[is.na(cor.hightoplot)] <- 0 
  cor.range <- c(quantile(as.vector(cor.hightoplot),probs=c(0,.5,.7,.85)),1)
  
  missing.summary <- missing.summary[match(names(High.cor),names(missing.summary))[1:min(numPlots,length(High.cor))]]
  
  ## nested functions -----
  
  ## panel histograms function ---------------------------------------------------
  
  # puts histograms on the diagonal
  
  panel.hist <- function(x, ...){
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5) )
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks; nB <- length(breaks)
    y <- h$counts; y <- y/max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col="steelblue", ...)
  }
  
  ## panel correlation function --------------------------------------------------
  
  # put (absolute) correlations on the upper panels,
  # with size proportional to the correlations.
  
  panel.cor <- function(x, y, 
                        digits=2, 
                        prefix="", 
                        cor.range,
                        cor.mult,
                        minCor,
                        ...){
    
    a <- colors()
    usr <- par("usr"); on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y,use="pairwise.complete.obs"))
    spear <- abs(cor(x,y,method="spearman",use="pairwise.complete.obs"))
    ken <- abs(cor(x,y,method="kendall",use="pairwise.complete.obs"))
    all.cor <- max(r,spear,ken)
    ramp <- heat.colors(20, alpha = 0.7)[20:1]
    if(all.cor >= minCor){
      rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col =
             ramp[which.min(abs(all.cor - seq(from=minCor,to=1,length=20)))])
    }
    r <- max(all.cor)
    cex.cor <- 3*cor.mult
    txt <- format(c(r, 0.123456789), digits=digits)[1]
    txt <- paste(prefix, txt, sep="")
    #if(missing(cex.cor)) cex.cor <- 1.2/strwidth(txt)
    
    txt2=""
    if(max(all.cor) > cor.range[2]){
      if(spear == max(all.cor) && spear != cor(x,y,use="pairwise.complete.obs")) {
        txt2 <- " s"
      } else if(ken == max(all.cor) && ken != cor(x,y,use="pairwise.complete.obs")){
        txt2 <-" k"
      }
      
    }
    text(0.5, 0.5, txt, cex = 0.7+cex.cor * (r-min(cor.range))/(max(cor.range)-min(cor.range)))
    text(0.9,0.1,txt2,cex=cor.mult)
  }
  
  # create output image -----
  
  options(warn=-1)
  numPlots <- min(ncol(HighToPlot),numPlots)
  if(numPlots<8){ 
    wdth = 1500
    cex.mult = 3
    } else if (numPlots<15){ wdth = 3000
        if(numPlots<12){ cex.mult=4 
        } else { cex.mult=3 }
      } else { 
        wdth = 4500
        if(numPlots<17){ cex.mult = 4
        } else { cex.mult=3 }
      }
  
  png(outputFile, width=wdth, height=wdth, pointsize=13)
 
  myPairs(for.dev = for.dev,
          minCor = minCor,
          missing.summary = missing.summary,
          my.labels = (as.vector(High.cor)[1:numPlots]),
          lower.panel = panel.smooth,
          diag.panel = panel.hist, 
          upper.panel = panel.cor,
          pch = 21,
          bg = c("blue","red","yellow")[factor(for.dev$response,levels=c(0,1,-9999))],
          col.smooth = "red",
          cex.mult = cex.mult,
          cor.range = cor.range,
          oma = c(0,2,6,0),
          family = family)
  
  graphics.off()
  options(warn=0)
  
}

## My Panel Smooth function ----------------------------------------------------

my.panel.smooth <- function(x,
                            y, 
                            col = par("col"),
                            bg = NA, 
                            pch = par("pch"),
                            family = "binomial",
                            cex = 1, 
                            col.smooth = "red",
                            span = 2/3, 
                            iter = 3, 
                            weights = rep(1,times=length(y)),
                            cex.mult,
                            Ylab,
                            plot.it = TRUE,
                            lin = 1,
                            ...)
{
  #This function fits a gam to show the relationship between a binary response and the specified predictor
  #similar to a lowess smooth but appropriate for binary response.  Occasionally gam fails or doesn't converge
  #this is indicated by the null deviance being less than the fit deviance.  When this occurs a glm is fit.

  o <- order(x)
  x <- x[o]
  y <- y[o]
  
  if(sum(y==0)/sum(y==1)>1.2 & family!="poisson"){ 
    wgt <- c(sum(y==1)/sum(y==0),1)[factor(y,levels=c(0,1))]
  } else { wgt <- rep(1,times=length(y)) }
  
  if(!is.factor(x)){
    options(warn=2)
    g <- try(gam(y~s(x,2),weights=wgt,family=family),silent=TRUE)
    options(warn=-1)
    dev.broke <- try((1-g$dev/g$null.deviance)<0,silent=TRUE)
    if(class(dev.broke)=="try-error"){ dev.broke=TRUE }
    if("try-error" %in% class(g) | dev.broke){
      gam.failed <- TRUE
      g <- glm(y~x+x^2,weights=wgt,family=family)
      y.fit <- predict(g, type="response")
    } else {
      y.fit <- predict.Gam(g, type="response")
      gam.failed <- FALSE
    }
  } else {
    g <- glm(y~x,weights=wgt,family=family)
    y.fit <- predict(g, type="response")
    if(plot.it){
      if(family == "poisson"){ 
        col.palatte <- c("blue",heat.colors(length(unique(y))))
        main.lab <- "Proportion of Count by Category"
      }
      else { 
        col.palatte <- c("blue","red")
        main.lab <- "Proportion of Pres/Abs by Category"
      }       
      barplot(prop.table(table(y,x),margin=2),main=main.lab,col=col.palatte,add=TRUE,...)
      mtext(paste("% dev exp ",round(100*(1-g$dev/g$null.deviance),digits=1),sep=""),side=2,line=lin,cex=cex.mult*.7)
    }
    return(100*(1-g$dev/g$null.deviance)) 
  } 
  if(plot.it){
    if(family == "poisson"){
      points(x, y, bg="black",col="black",...)
    } else {
      points(x, y, pch = pch,
             bg = c("blue","red")[factor(y,levels=c(0,1))],
             col = c("blue4","red4")[factor(y,levels=c(0,1))],
             cex=cex.mult,...)
    } 
    segments(x0=x[1:(length(x)-1)],
             y0=y.fit[1:(length(x)-1)],
             x1=x[2:length(x)],
             y1=y.fit[2:length(x)],
             col="red",
             cex=3*cex.mult,lwd=cex.mult)
    if(missing(Ylab)){
      mtext(paste("% dev exp ",round(100*(1-g$dev/g$null.deviance),digits=1),sep=""),side=2,line=lin,cex=cex.mult*.7)
    }
    return(gam.failed)
  } else { return(100*(1-g$dev/g$null.deviance)) }
  
}


## My Pairs function -----------------------------------------------------------

myPairs <- function(for.dev,
                    minCor,
                    missing.summary,
                    my.labels,
                    labels = NULL, 
                    panel = points,
                    lower.panel, # = panel,
                    upper.panel, # = panel,
                    diag.panel, # = NULL, 
                    text.panel = textPanel,
                    label.pos = 0.5 + has.diag/3, 
                    cex.labels = NULL, 
                    font.labels = 1,
                    row1attop = TRUE, 
                    gap = 1,
                    cex.mult,
                    family="binomial",
                    oma = NULL,
                    main = NULL, 
                    ...)
{
  
  ### Nested functions -----
  
  textPanel <- function(x = 0.5, y = 0.5, txt, cex, font){
    text(x, y, txt, cex = cex, font = font)
  }
  
  localAxis <- function(side, x, y, xpd, bg, col = NULL, main, oma, ...) {
    if (side %% 2 == 1){
      Axis(x, side = side, xpd = NA, ...)
    } else {
      Axis(y, side = side, xpd = NA, ...)
    }
  }
  
  localPlot <- function(..., main, oma, font.main, cex.main) { plot(...) }
  localLowerPanel <- function(..., main, oma, font.main, cex.main){lower.panel(...)} 
  localUpperPanel <-  function(..., main, oma, font.main, cex.main){upper.panel(...)}
  localDiagPanel <- function(..., main, oma, font.main, cex.main){ diag.panel(...) }
  
  # code starts for my pairs function ----
  
  response <- for.dev$response
  x <- for.dev$dat
  
  if (!is.matrix(x)) {
    x <- as.data.frame(x)
    for (i in seq_along(names(x))) {
      if (is.factor(x[[i]]) || is.logical(x[[i]])){ x[[i]] <- as.numeric(x[[i]]) }
      if (!is.numeric(unclass(x[[i]]))){ stop("non-numeric argument to 'pairs'") }
    }
  } else if (!is.numeric(x)) { stop("non-numeric argument to 'pairs'") }
    
  panel <- match.fun(panel)
  if ((has.lower <- !is.null(lower.panel)) && !missing(lower.panel)){ lower.panel <- match.fun(lower.panel) }
  if ((has.upper <- !is.null(upper.panel)) && !missing(upper.panel)){ upper.panel <- match.fun(upper.panel) }
  if ((has.diag <- !is.null(diag.panel)) && !missing(diag.panel)){ diag.panel <- match.fun(diag.panel) }

  if (row1attop) {
    tmp <- lower.panel
    lower.panel <- upper.panel
    upper.panel <- tmp
    tmp <- has.lower
    has.lower <- has.upper
    has.upper <- tmp
  }
  
  nc <- ncol(x)
  if (nc < 2){ stop("only one column in the argument to 'pairs'") }
   
  has.labs <- TRUE
  labels <- colnames(x)

  if (is.null(oma)) {
    oma <- c(4, 4, 4, 4)
    if (!is.null(main))
      oma[3L] <- 6
  }
  nCol <- ifelse(length(unique(response))>1, nc+1, nc)
  j.start <- ifelse(length(unique(response))>1, 0, 1)
  opar <- par(mfrow = c(nc, nCol), mar = rep.int(gap/2, 4))
  on.exit(par(opar))
  par(oma=oma)
  for (i in if (row1attop){ 1L:(nc) } else {nc:1L}) 
    for (j in j.start:(nc)) {
      
      top.gap <- c(rep(gap/2,times=nc))
      bottom.gap <- c(rep(gap/2,times=nc-1),3*gap)
      left.gap <- c(3*gap,3*gap,rep(gap/2,times=nc-1))
      par(mar = c(bottom.gap[i],left.gap[j+1],top.gap[i],gap/2))
      
      if(j==0){
        localPlot(x[, i], response, xlab = "", ylab = "", axes = FALSE, type="n",...)
        if(i==1){ mtext("Response",line=.3,cex=.7*cex.mult) }
        box()
        my.lab <- paste("cor=",round(max(abs(cor(x[,(i)],response,use="pairwise.complete.obs")),
                                         abs(cor(x[,(i)],response,method="spearman",use="pairwise.complete.obs")),
                                         abs(cor(x[,(i)],response,method="kendall",use="pairwise.complete.obs"))),digits=2),sep="")
        if(family == "gaussian"){ 
          panel.smooth(as.vector(x[, (i)]), as.vector(response),...)
          title(ylab=paste("cor=",round(max(abs(cor(x[,(i)],response,use="pairwise.complete.obs")),
                                            abs(cor(x[,(i)],response,method="spearman",use="pairwise.complete.obs")),
                                            abs(cor(x[,(i)],response,method="kendall",use="pairwise.complete.obs"))),digits=2),
                           sep=""),line=.02,cex.lab=1.5)
        } else if(missing(for.dev)){
          pct.dev <- try(my.panel.smooth(x = as.vector(x[, (i)]), 
                                         y = response,
                                         cex.mult=cex.mult,
                                         cex.lab=cex.mult,
                                         line=1,
                                         family=family,...),silent=TRUE)
        } else {
          pct.dev <- try(my.panel.smooth(as.vector(for.dev$dat[, (i)]),
                                         as.vector(for.dev$response),
                                         cex.mult=cex.mult,
                                         cex.lab=cex.mult,
                                         line=1,
                                         family=family,...),silent=TRUE)
        }     
      } else{
        
        localPlot(x[, j], x[, i], xlab = "", ylab = "", axes = FALSE, type="n",...)
        
        if (i == j || (i < j && has.lower) || (i > j && has.upper)) {
          box()
          if(i==1) {
            mtext(paste("Total Cor=",my.labels[j],sep=""),side=3,line=ifelse(missing.summary[j]>0.03,3,0.3),cex=0.65*cex.mult)
            if(missing.summary[j]>0.03){ mtext(paste(round(missing.summary[j]*100), "% missing",sep=""),side=3,line=0.3,cex=cex.mult*0.55) }
          }
          if (i == nc) { localAxis(3 - 2 * row1attop, x[, j], x[, i],cex.axis=cex.mult*0.5, ...) }
          if (j == 1 && (i!=1 || !has.upper || !has.lower)) { localAxis(2, x[, j], x[, i],cex.axis=cex.mult*0.5, ...) }
          
          mfg <- par("mfg")
          if (i == j) {
            if (has.diag)
              localDiagPanel(as.vector(x[, i]),...)
            if (has.labs) {
              par(usr = c(0, 1, 0, 1))
              if(i==1){
                for(k in 1:length(labels)){
                  if((lng<-nchar(labels[k]))>=14) labels[k]<-paste(substr(labels[k],1,12),"\n",substr(labels[k],13,lng),sep="")
                }
                if (is.null(cex.labels)) {
                  l.wid <- strwidth(labels, "user")
                  cex.labels <- max(0.8, min(2, 0.9/max(l.wid)))
                }
              }
              text.panel(0.5, label.pos, labels[i], cex = 0.45*cex.labels*cex.mult,
                         font = font.labels)
            }
          }
          else if (i < j)
            if(length(unique(x[,i])>2)){
              localLowerPanel(as.vector(x[, j]), as.vector(x[,i]),cex=cex.mult*3,cor.mult=cex.mult,minCor=minCor,...) 
              } else {
                if(missing(for.dev)){ 
                  pct.dev <- try(my.panel.smooth(as.vector(x[,j]), 
                                                 as.vector(x[,i]),
                                                 cex.mult=cex.mult,
                                                 cex.lab=cex.mult,
                                                 line=1,
                                                 family=family,...),silent=TRUE)
                } else {
                  pct.dev <- try(my.panel.smooth(as.vector(for.dev$dat[, (i)]),
                                                 as.vector(for.dev$response),
                                                 cex.mult=cex.mult,
                                                 cex.lab=cex.mult,
                                                 line=1,
                                                 family=family,...),silent=TRUE)
                  }} else {
                    localUpperPanel(as.vector(x[, j]), 
                                    as.vector(x[,i]),
                                    cex=cex.mult,
                                    cor.mult=cex.mult, 
                                    minCor = minCor,
                                    ...)
                    }
          if (any(par("mfg") != mfg)){ stop("the 'panel' function made a new plot") }
        } else { par(new = FALSE) }
      }}
  if (!is.null(main)) {
    if (is.null(font.main)){ font.main <- par("font.main") }
    if (is.null(cex.main)){ cex.main <- par("cex.main") }
    mtext(main, 3, 3, TRUE, 0.5, cex = cex.main, font = font.main)
  }
  invisible(NULL)
}


