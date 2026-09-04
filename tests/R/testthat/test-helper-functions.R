## Unit tests for src/00-helper-functions.R -- pure prediction, evaluation, and
## data-wrangling helpers that don't require a live SyncroSim scenario.

test_that("calc.deviance computes binomial deviance correctly", {
  source_wisdm_r("00-helper-functions.R")

  obs <- c(1, 0, 1, 0)
  preds <- c(0.9, 0.1, 0.8, 0.2)
  expected <- -2 * mean(obs * log(preds) + (1 - obs) * log(1 - preds))

  expect_equal(calc.deviance(obs, preds, family = "binomial"), expected)
})

test_that("calc.deviance computes gaussian deviance correctly", {
  source_wisdm_r("00-helper-functions.R")

  obs <- c(1, 2, 3)
  preds <- c(1.5, 2.5, 2.5)

  expect_equal(
    calc.deviance(obs, preds, family = "gaussian", calc.mean = FALSE),
    sum((obs - preds)^2)
  )
})

test_that("calc.deviance errors on mismatched lengths", {
  source_wisdm_r("00-helper-functions.R")
  expect_error(calc.deviance(c(1, 0), c(0.5)), "equal length")
})

test_that("glm.predict and pred.fct agree and return response-scale predictions", {
  source_wisdm_r("00-helper-functions.R")

  set.seed(1)
  x1 <- rnorm(50)
  y <- rbinom(50, 1, plogis(0.5 * x1))
  x <- data.frame(x1 = x1)
  model <- glm(y ~ x1, data = data.frame(y = y, x1 = x1), family = binomial)

  direct <- glm.predict(model, x)
  viaPredFct <- pred.fct(model, x, modType = "glm")

  expect_equal(direct, viaPredFct)
  expect_true(all(direct >= 0 & direct <= 1))
})

test_that("safe_rbind aligns mismatched columns and fills NA for the rest", {
  source_wisdm_r("00-helper-functions.R")

  df <- data.frame(a = 1:2, b = c("x", "y"), stringsAsFactors = FALSE)
  newRow <- data.frame(b = "z", c = 10, stringsAsFactors = FALSE)

  result <- safe_rbind(df, newRow)

  expect_equal(nrow(result), 3)
  expect_equal(names(result), c("a", "b", "c"))
  expect_true(is.na(result$a[3]))
  expect_true(is.na(result$c[1]))
  expect_equal(result$c[3], 10)
})

test_that("calcSiteWeights weights background points by the presence:absence ratio", {
  source_wisdm_r("00-helper-functions.R")

  response <- c(1, 1, 1, 0, 0, 0, 0, 0, 0) # 3 presence, 6 absence
  weights <- calcSiteWeights(response)

  expect_equal(weights[response == 1], rep(1, 3))
  expect_equal(unique(weights[response == 0]), 3 / 6)
})

test_that("boyceIndex is high for well-calibrated predictions and NA for degenerate input", {
  source_wisdm_r("00-helper-functions.R")

  set.seed(42)
  allPreds <- runif(500)
  # presence predictions skewed toward the high end of the suitability range
  presPreds <- pmin(allPreds[1:200] + 0.4, 1)

  bi <- boyceIndex(allPreds, presPreds)
  expect_true(is.numeric(bi))
  expect_true(bi > 0)

  expect_true(is.na(boyceIndex(runif(2), runif(2))))
})

test_that("modelEvaluation reports AUC of 1 for perfectly separated predictions", {
  source_wisdm_r("00-helper-functions.R")

  predOcc <- c(0.9, 0.8, 0.7)
  predAbs <- c(0.3, 0.2, 0.1)

  ev <- modelEvaluation(predOcc, predAbs)

  expect_s4_class(ev, "modelEvaluation")
  expect_equal(ev@stats$auc, 1)
  expect_equal(ev@stats$np, 3L)
  expect_equal(ev@stats$na, 3L)
})

test_that("modelEvaluation errors when presence or absence data is empty", {
  source_wisdm_r("00-helper-functions.R")
  expect_error(modelEvaluation(numeric(0), c(0.1, 0.2)), "cannot evaluate")
})
