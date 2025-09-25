## -------------------------
## wisdm - helper functions
## ApexRMS, March 2022
## -------------------------

# Calculate Deviance function --------------------------------------------------

calc.deviance <- function(
  obs, # observed response
  preds, # predicted response
  weights = rep(1, length(obs)),
  family = "binomial",
  calc.mean = TRUE,
  return.list = FALSE
) {
  # function to calculate deviance given two vectors of raw and fitted values
  # requires a family argument which is set to binomial by default

  if (length(obs) != length(preds)) {
    stop("observations and predictions must be of equal length")
  }

  if (family == "binomial" | family == "bernoulli") {
    # preds[preds==0] <- 0.000000001
    deviance.contribs <- (obs * log(preds)) + ((1 - obs) * log(1 - preds))
    deviance <- -2 * sum(deviance.contribs * weights, na.rm = T)
  }

  if (family == "laplace") {
    deviance <- sum(abs(obs - preds))
  }

  if (family == "gaussian") {
    deviance <- sum((obs - preds) * (obs - preds))
  }

  if (calc.mean) {
    deviance <- deviance / length(obs)
  }

  if (return.list) {
    deviance <- list(deviance = deviance, dev.cont = deviance.contribs)
  }
  return(deviance)
}

# Predict Function -------------------------------------------------------------

pred.fct <- function(
  mod, # mod = the model fit object
  x, # x = data to predict for
  modType
) {
  # modType = one of mars, glm, rf, brt, maxlike at present

  y <- rep(NA, nrow(x))
  idx <- stats::complete.cases(x)

  if (modType == "glm") {
    y <- tryCatch(glm.predict(mod, x), error = function(e) {
      rep(NA, nrow(x))
    })
  }
  if (modType == "rf") {
    preds <- tryCatch(
      rf.predict(mod, x[idx, , drop = FALSE]),
      error = function(e) {
        rep(NA, sum(idx))
      }
    )
    y[idx] <- preds
  }
  if (modType == "maxent") {
    preds <- tryCatch(
      maxent.predict(mod, x[idx, , drop = FALSE]),
      error = function(e) {
        rep(NA, sum(idx))
      }
    )
    y[idx] <- preds
  }
  if (modType == "brt") {
    preds <- tryCatch(
      brt.predict(mod, x[idx, , drop = FALSE]),
      error = function(e) {
        rep(NA, sum(idx))
      }
    )
    y[idx] <- preds
  }
  if (modType == "gam") {
    preds <- tryCatch(
      gam.predict(mod, x[idx, , drop = FALSE]),
      error = function(e) {
        rep(NA, sum(idx))
      }
    )
    y[idx] <- preds
  }
  return(y)
} # end pred.vals function

## glm predict function --------------------------------------------------------

glm.predict <- function(model, x) {
  # retrieve key items from the global environment #
  # make predictions.

  y <- as.vector(stats::predict(object = model, newdata = x, type = "response"))

  # encode missing values as -1.
  y[is.na(y)] <- NaN

  # return predictions.
  return(y)
}

## rf predict function ---------------------------------------------------------

rf.predict <- function(model, x) {
  # retrieve key items from the global environment #
  # make predictions from complete data only #
  y <- rep(NA, nrow(x))
  idx <- stats::complete.cases(x)
  if (any(idx)) {
    y[idx] <- tryCatch(
      as.vector(predict(
        model,
        newdata = x[idx, , drop = FALSE],
        type = "vote"
      )[, 2]),
      error = function(e) {
        rep(NA, sum(idx))
      }
    )
  }

  # return predictions.
  return(y)
}

## tweak prediction function (for rf) ------------------------------------------

tweak.p <- function(p) {
  p[p == 1] <- max(p[p < 1])
  p[p == 0] <- min(p[p > 0])
  return(p)
}

## maxent predict function -----------------------------------------------------

maxent.predict <- function(model, x) {
  if (names(model[[1]])[1] == "Raw.coef") {
    attach(model[[1]])
    on.exit(detach(model[[1]]))
  } else {
    attach(model)
    on.exit(detach(model))
  }

  # build prediction shell filled with NA
  prediction <- rep(NA, nrow(x))
  names(prediction) <- row.names(x)

  # These all default to zero in case the feature was excluded
  Raw <- Quad <- Prod <- Forw <- Rev <- Thresh <- 0

  if (!is.null(Raw.coef)) {
    Raw <- model.matrix(
      as.formula(paste("~", paste(Raw.coef[, 1], collapse = " + "), sep = "")),
      x
    )
    Raw <- apply(t(Raw) * Raw.mult, 2, sum)
  }
  if (!is.null(Quad.coef)) {
    Quad <- model.matrix(
      as.formula(paste(
        "~ I(",
        paste(Quad.coef[, 1], collapse = ") + I("),
        ")",
        sep = ""
      )),
      x
    )
    Quad <- apply(t(Quad) * Quad.mult, 2, sum)
  }
  if (!is.null(Prod.coef)) {
    Prod <- model.matrix(
      as.formula(paste(
        "~ ",
        paste(Prod.coef[, 1], collapse = " + "),
        sep = ""
      )),
      x
    )
    Prod <- apply(t(Prod) * Prod.mult, 2, sum)
  }
  if (!is.null(Fwd.Hinge)) {
    Forw <- model.matrix(
      as.formula(paste(
        "~ ",
        paste(
          "-1 + ",
          paste(
            "I(",
            paste(
              Fwd.Hinge[, 1],
              paste(
                Fwd.Hinge[, 1],
                paste("(", Fwd.Hinge[, 3], ")", sep = ""),
                sep = ">"
              ),
              sep = "*("
            ),
            collapse = ")) + "
          ),
          "))",
          sep = ""
        ),
        sep = ""
      )),
      x
    )
    Forw.Cst <- Forw != 0
    Forw <- apply(t(Forw) * FH.mult, 2, sum) +
      apply(t(Forw.Cst) * FH.cnst, 2, sum)
  }
  if (!is.null(Rev.Hinge)) {
    Rev <- model.matrix(
      as.formula(paste(
        "~ ",
        paste(
          "-1 + ",
          paste(
            "I(",
            paste(
              Rev.Hinge[, 1],
              paste(
                Rev.Hinge[, 1],
                paste("(", Rev.Hinge[, 4], ")", sep = ""),
                sep = "<"
              ),
              sep = "*("
            ),
            collapse = ")) + "
          ),
          "))",
          sep = ""
        ),
        sep = ""
      )),
      x
    )
    Rev.Cst <- Rev != 0
    Rev <- apply(t(Rev) * Rev.mult, 2, sum) +
      apply(t(Rev.Cst) * Rev.cnst, 2, sum)
  }
  if (!is.null(Thresh.val)) {
    Thresh.val[, 1] <- gsub("=", "==", Thresh.val[, 1])
    Thresh <- model.matrix(
      as.formula(paste(
        "~ ",
        paste("I", Thresh.val[, 1], collapse = " + ", sep = ""),
        sep = ""
      )),
      x
    )
    Thresh <- apply(t(Thresh[, 2:ncol(Thresh)]) * Thresh.cnst, 2, sum)
  }

  S <- Raw + Quad + Prod + Forw + Rev + Thresh - normalizers[1, 2]

  qx <- exp(S) / normalizers[2, 2]
  predVals <- qx * exp(entropy) / (1 + qx * exp(entropy))
  prediction[names(predVals)] <- predVals
  return(prediction)
}

## brt predict function --------------------------------------------------------

brt.predict <- function(model, x) {
  y <- rep(NA, nrow(x))
  idx <- stats::complete.cases(x)
  if (any(idx)) {
    ntrees <- if (!is.null(model$gbm.call$best.trees)) {
      model$gbm.call$best.trees
    } else if (!is.null(model$n.trees)) {
      model$n.trees
    } else {
      stop("Number of trees not found in BRT model")
    }
    y[idx] <- gbm::predict.gbm(
      model,
      x[idx, , drop = FALSE],
      n.trees = ntrees,
      type = "response"
    )
  }

  # # make predictions from full data #
  # y <- predict.gbm(model, x, model$target.trees, type="response")
  # # encode missing values as -1.
  # a <- complete.cases(x)
  # y[!(a)]<- NaN

  # return predictions.
  return(y)
}

## gam predict function --------------------------------------------------------

gam.predict <- function(model, x) {
  y <- rep(NA, nrow(x))
  idx <- stats::complete.cases(x)
  y[idx] <- mgcv::predict.gam(model, x[idx, ], type = "response")

  # return predictions.
  return(y)
}

# Evaluate Function ------------------------------------------------------------

#' Model Evaluation Function
#'
#' This function takes in the predicted values from a species distribution model
#' (predicted presences and absences) and returns an object of class 'modelEvaluation'
#' that contains various evaluation metrics such as the confusion matrix and threshold values.
#'
#' @param predOcc A numeric vector of predicted values for observed presences from a species distribution model.
#' @param predAbs A numeric vector of predicted values for observed absences from a species distribution model.
#'
#' @return An object of class 'modelEvaluation' with the following slots:
#'   * presence: The input presence values after removing NAs.
#'   * absence: The input absence values after removing NAs.
#'   * confusion: A confusion matrix with columns for true positives (tp), false positives (fp),
#'     false negatives (fn), and true negatives (tn) for each threshold.
#'   * stats: A data frame with overall evaluation metrics including the number of presence and
#'     absence points, prevalence, AUC, correlation, p-value of the correlation, and overall diagnostic power.
#'   * tr_stats: A data frame with threshold-specific evaluation metrics including kappa, correct classification rate (CCR),
#'     true positive rate (TPR), true negative rate (TNR), false positive rate (FPR), false negative rate (FNR),
#'     positive predictive power (PPP), negative predictive power (NPP), misclassification rate (MCR), and odds ratio (OR).
#'   * thresholds: A data frame with various threshold values including the threshold with maximum kappa,
#'     the threshold with maximum sum of sensitivity and specificity, the threshold with no omission,
#'     the threshold with equal prevalence, and the threshold with equal sensitivity and specificity.
#'
#' @examples
#' # Assuming predOcc and predAbs are your predicted presence and absence values
#' result <- modelEvaluation(predOcc, predAbs)

modelEvaluation <- function(predOcc, predAbs) {
  p <- stats::na.omit(predOcc)
  a <- stats::na.omit(predAbs)
  np <- length(p)
  na <- length(a)
  if (na == 0 | np == 0) {
    stop(
      "cannot evaluate a model without absence and presence data that are not NA"
    )
  }

  if (TRUE) {
    if (length(p) > 1000) {
      tr <- as.vector(stats::quantile(p, 0:1000 / 1000))
    } else {
      tr <- p
    }
    if (length(a) > 1000) {
      tr <- c(tr, as.vector(stats::quantile(a, 0:1000 / 1000)))
    } else {
      tr <- c(tr, a)
    }
    tr <- sort(unique(round(tr, 8)))
    tr <- c(tr - 0.0001, tr[length(tr)] + c(0, 0.0001))
  } else {
    tr <- sort(as.vector(tr))
  }

  N <- na + np

  if (!methods::isClass("modelEvaluation")) {
    setClass(
      "modelEvaluation",
      slots = list(
        presence = "numeric",
        absence = "numeric",
        confusion = "matrix",
        stats = "data.frame",
        tr_stats = "data.frame",
        thresholds = "data.frame"
      )
    )
  }

  xc <- methods::new("modelEvaluation")
  xc@presence = p
  xc@absence = a

  R <- sum(rank(c(p, a))[1:np]) - (np * (np + 1) / 2)
  auc <- R / (as.numeric(na) * as.numeric(np))

  cr <- try(
    stats::cor.test(c(p, a), c(rep(1, length(p)), rep(0, length(a)))),
    silent = TRUE
  )
  corc <- pcor <- NA
  if (!inherits(cr, "try-error")) {
    corc <- cr$estimate
    pcor <- cr$p.value
  }

  res <- matrix(ncol = 4, nrow = length(tr))
  colnames(res) <- c("tp", "fp", "fn", "tn")
  for (i in 1:length(tr)) {
    res[i, 1] <- length(p[p >= tr[i]]) # a  true positives
    res[i, 2] <- length(a[a >= tr[i]]) # b  false positives
    res[i, 3] <- length(p[p < tr[i]]) # c  false negatives
    res[i, 4] <- length(a[a < tr[i]]) # d  true negatives
  }
  xc@confusion = res
  a = res[, 1]
  b = res[, 2]
  c = res[, 3]
  d = res[, 4]
  # after Fielding and Bell
  np <- as.integer(np)
  na <- as.integer(na)
  prevalence = (a[1] + c[1]) / N
  # overall diagnostic power
  ODP = (b[1] + d[1]) / N
  xc@stats <- data.frame(np, na, prevalence, auc, cor = corc, pcor, ODP)
  rownames(xc@stats) <- NULL
  # correct classification rate
  CCR = (a + d) / N
  # sensitivity, or true positive rate
  TPR = a / (a + c)
  # specificity, or true negative rate
  TNR = d / (b + d)
  # False positive rate
  FPR = b / (b + d)
  # False negative rate
  FNR = c / (a + c)
  PPP = a / (a + b)
  NPP = d / (c + d)
  # misclassification rate
  MCR = (b + c) / N
  # odds ratio
  OR = (a * d) / (c * b)

  prA = (a + d) / N
  prY = (a + b) / N * (a + c) / N
  prN = (c + d) / N * (b + d) / N
  prE = prY + prN
  kappa = (prA - prE) / (1 - prE)
  xc@tr_stats <- data.frame(
    treshold = tr,
    kappa,
    CCR,
    TPR,
    TNR,
    FPR,
    FNR,
    PPP,
    NPP,
    MCR,
    OR
  )

  max_kappa <- tr[which.max(kappa)]
  # maximum sum of the sensitivity (true positive rate) and specificity (true negative rate)
  max_spec_sens <- tr[which.max(TPR + TNR)]
  # no omission
  no_omission <- tr[max(which(res[, "fn"] == 0))]
  # Suggestions by Diego Nieto-Lugilde
  # equal prevalence
  equal_prevalence = tr[which.min(abs(tr - prevalence))]
  # equal sensitivity and specificity
  equal_sens_spec <- tr[which.min(abs(TPR - TNR))]

  xc@thresholds <- data.frame(
    max_kappa,
    max_spec_sens,
    no_omission,
    equal_prevalence,
    equal_sens_spec
  )

  return(xc)
}

# safe rbind function -------------------------------------------------------

safe_rbind <- function(df, row) {
  # Ensure row is a data.frame
  if (!is.data.frame(row)) {
    row <- as.data.frame(row, stringsAsFactors = FALSE)
  }

  # Add any missing columns to row
  for (col in setdiff(names(df), names(row))) {
    row[[col]] <- NA
  }

  # Add any new columns to df
  for (col in setdiff(names(row), names(df))) {
    df[[col]] <- NA
  }

  # Align column order
  row <- row[names(df)]

  # Bind
  df <- rbind(df, row)
  rownames(df) <- NULL # reset row names
  return(df)
}
