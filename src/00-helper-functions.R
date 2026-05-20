## -------------------------
## wisdm - helper functions
## ApexRMS, March 2022
## -------------------------

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "00-constants.R"))

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

  predictSafe <- function(predictFct, mod, x, idx) {
    tryCatch(predictFct(mod, x[idx, , drop = FALSE]), error = function(e) {
      rep(NA, sum(idx))
    })
  }

  if (modType == "glm") {
    y[idx] <- predictSafe(glm.predict, mod, x, idx)
  }
  if (modType == "rf") {
    y[idx] <- predictSafe(rf.predict, mod, x, idx)
  }
  if (modType == "maxent") {
    y[idx] <- predictSafe(maxent.predict, mod, x, idx)
  }
  if (modType == "brt") {
    y[idx] <- predictSafe(brt.predict, mod, x, idx)
  }
  if (modType == "gam") {
    y[idx] <- predictSafe(gam.predict, mod, x, idx)
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
  ## Helper to build model matrices safely (no intercepts)
  makeMatrix <- function(formula_txt, data) {
    form <- as.formula(formula_txt)
    mm <- model.matrix(
      form,
      model.frame(form, data = data, na.action = na.pass)
    )
    mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
    mm
  }

  ## Helper for feature blocks with possible intercept
  featurePredict <- function(mat, mult, name) {
    if (is.null(mult) || is.null(mat)) {
      return(0)
    }
    if (length(mult) == ncol(mat) + 1) {
      warning(paste0(
        name,
        ": assuming mult[1] is intercept; verify model structure"
      ))
      as.vector(mat %*% mult[-1]) + mult[1]
    } else if (length(mult) == ncol(mat)) {
      as.vector(mat %*% mult)
    } else {
      stop(paste0(name, ".mult length mismatch with ", name, ".mat columns"))
    }
  }

  # Attach the model structure (handles list-wrapped models)
  if (names(model[[1]])[1] == "Raw.coef") {
    m <- model[[1]]
  } else {
    m <- model
  }

  # Initialize components
  n <- nrow(x)
  Raw <- Quad <- Prod <- Forw <- Rev <- Thresh <- rep(0, n)

  # Raw Features ---
  if (!is.null(m$Raw.coef)) {
    Raw.mat <- makeMatrix(
      paste("~ -1 + ", paste(m$Raw.coef[, 1], collapse = " + ")),
      x
    )
    Raw <- featurePredict(Raw.mat, m$Raw.mult, "Raw")
  }

  # Quadratic Features ---
  if (!is.null(m$Quad.coef)) {
    Quad.mat <- makeMatrix(
      paste0("~ -1 + ", paste0("I(", m$Quad.coef[, 1], ")", collapse = " + ")),
      x
    )
    Quad <- featurePredict(Quad.mat, m$Quad.mult, "Quad")
  }

  # Product Features ---
  if (!is.null(m$Prod.coef)) {
    Prod.mat <- makeMatrix(
      paste("~ -1 + ", paste(m$Prod.coef[, 1], collapse = " + ")),
      x
    )
    Prod <- featurePredict(Prod.mat, m$Prod.mult, "Prod")
  }

  # Forward Hinge ---
  if (!is.null(m$Fwd.Hinge)) {
    Forw.mat <- makeMatrix(
      paste0(
        "~ -1 + ",
        paste0(
          "I(",
          m$Fwd.Hinge[, 1],
          " * (",
          m$Fwd.Hinge[, 1],
          " > (",
          m$Fwd.Hinge[, 3],
          ")))",
          collapse = " + "
        )
      ),
      x
    )
    Forw.ind <- makeMatrix(
      paste0(
        "~ -1 + ",
        paste0(
          "as.numeric(",
          m$Fwd.Hinge[, 1],
          " > (",
          m$Fwd.Hinge[, 3],
          "))",
          collapse = " + "
        )
      ),
      x
    )

    if (ncol(Forw.mat) != length(m$FH.mult)) {
      stop("FH.mult length mismatch with Forw.mat columns")
    }
    if (ncol(Forw.ind) != length(m$FH.cnst)) {
      stop("FH.cnst length mismatch with Forw.ind columns")
    }

    Forw <- as.vector(Forw.mat %*% m$FH.mult) +
      as.vector(Forw.ind %*% m$FH.cnst)
  }

  # Reverse Hinge ---
  if (!is.null(m$Rev.Hinge)) {
    Rev.mat <- makeMatrix(
      paste0(
        "~ -1 + ",
        paste0(
          "I(",
          m$Rev.Hinge[, 1],
          " * (",
          m$Rev.Hinge[, 1],
          " < (",
          m$Rev.Hinge[, 4],
          ")))",
          collapse = " + "
        )
      ),
      x
    )
    Rev.ind <- makeMatrix(
      paste0(
        "~ -1 + ",
        paste0(
          "as.numeric(",
          m$Rev.Hinge[, 1],
          " < (",
          m$Rev.Hinge[, 4],
          "))",
          collapse = " + "
        )
      ),
      x
    )

    if (ncol(Rev.mat) != length(m$Rev.mult)) {
      stop("Rev.mult length mismatch with Rev.mat columns")
    }
    if (ncol(Rev.ind) != length(m$Rev.cnst)) {
      stop("Rev.cnst length mismatch with Rev.ind columns")
    }

    Rev <- as.vector(Rev.mat %*% m$Rev.mult) + as.vector(Rev.ind %*% m$Rev.cnst)
  }

  # Threshold Features ---
  if (!is.null(m$Thresh.val)) {
    thr_expr <- gsub(
      "(?<!(<|>|!|=))=(?!=)",
      "==",
      m$Thresh.val[, 1],
      perl = TRUE
    )
    Thresh.mat <- makeMatrix(
      paste("~ -1 + ", paste0("as.numeric(", thr_expr, ")", collapse = " + ")),
      x
    )
    if (ncol(Thresh.mat) != length(m$Thresh.cnst)) {
      stop("Thresh.cnst length mismatch with Thresh.mat columns")
    }
    Thresh <- as.vector(Thresh.mat %*% m$Thresh.cnst)
  }

  # Combine all feature groups
  S <- Raw + Quad + Prod + Forw + Rev + Thresh - m$normalizers[1, 2]

  # Apply logistic transform
  qx <- exp(S) / m$normalizers[2, 2]
  predVals <- qx * exp(m$entropy) / (1 + qx * exp(m$entropy))

  return(predVals)
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

# Launch a Shiny app with browser detection ------------------------------------

launchShinyApp <- function(appPath, appName = "Viewer") {
  # Search PATH first — works on all platforms (Linux, macOS, Windows)
  chrome.names <- c(
    "google-chrome",
    "google-chrome-stable",
    "chromium-browser",
    "chromium",
    "brave-browser",
    "microsoft-edge",
    "msedge"
  )
  firefox.names <- c("firefox")

  find.in.path <- function(candidates) {
    for (candidate in candidates) {
      found <- Sys.which(candidate)
      if (nzchar(found)) return(found)
    }
    NULL
  }

  browser.path <- find.in.path(chrome.names)
  browser.type <- if (!is.null(browser.path)) "chrome" else NULL

  if (is.null(browser.path)) {
    browser.path <- find.in.path(firefox.names)
    browser.type <- if (!is.null(browser.path)) "firefox" else NULL
  }

  # Windows fallback: browsers are often not on PATH, check common install dirs
  if (is.null(browser.path) && .Platform$OS.type == "windows") {
    local.app <- Sys.getenv("LOCALAPPDATA")
    prog.files <- c(Sys.getenv("ProgramFiles"), Sys.getenv("ProgramFiles(x86)"))
    chrome.paths <- c(
      file.path(local.app, "Google/Chrome/Application/chrome.exe"),
      file.path(prog.files, "Google/Chrome/Application/chrome.exe"),
      file.path(
        prog.files,
        "BraveSoftware/Brave-Browser/Application/brave.exe"
      ),
      file.path(prog.files, "Microsoft/Edge/Application/msedge.exe")
    )
    firefox.paths <- file.path(prog.files, "Mozilla Firefox/firefox.exe")

    for (path in chrome.paths) {
      if (file.exists(path)) {
        browser.path <- path
        browser.type <- "chrome"
        break
      }
    }
    if (is.null(browser.path)) {
      for (path in firefox.paths) {
        if (file.exists(path)) {
          browser.path <- path
          browser.type <- "firefox"
          break
        }
      }
    }
  }

  # On headless Unix (Docker or bare Linux server with no DISPLAY), bind to
  # 0.0.0.0:3838 so the app is reachable via port forwarding from the host.
  # A browser is still launched if one is found (e.g. via X11 forwarding).
  # Chrome/Firefox require --no-sandbox when running as root.
  # On non-headless Unix or Windows, behaviour is unchanged.
  headless_unix <- .Platform$OS.type == "unix" &&
    nchar(Sys.getenv("DISPLAY")) == 0

  if (headless_unix) {
    port <- 3838
    progressBar(
      type = "message",
      message = paste0(
        "ACTION REQUIRED: Open http://localhost:",
        port,
        " in your browser to use the Interactive ",
        appName
      )
    )
    shiny::runApp(
      appDir = appPath,
      host = "0.0.0.0",
      port = port,
      launch.browser = FALSE
    )
  } else if (is.null(browser.path)) {
    shiny::runApp(appDir = appPath, launch.browser = TRUE)
  } else {
    shiny::runApp(appDir = appPath, launch.browser = function(shinyurl) {
      if (browser.type == "chrome") {
        system2(
          browser.path,
          args = c(paste0("--app=", shinyurl), "--incognito"),
          wait = FALSE
        )
      } else {
        system2(
          browser.path,
          args = c("--private-window", shinyurl),
          wait = FALSE
        )
      }
    })
  }
}

# Calculate presence-weighted case weights -------------------------------------

calcSiteWeights <- function(response) {
  prNum <- as.numeric(table(response)["1"])
  bgNum <- as.numeric(table(response)["0"])
  ifelse(response == 1, 1, prNum / bgNum)
}
