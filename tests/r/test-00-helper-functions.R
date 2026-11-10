## -------------------------
## wisdm - tests for src/00-helper-functions.R
## ApexRMS, October 2026
## -------------------------

# calc.deviance ----------------------------------------------------------------

test_that("calc.deviance matches a hand-calculated binomial deviance", {
  siteData <- read.csv(test_path("fixtures", "deviance-example.csv"))

  # Presences contribute log(p), absences log(1 - p):
  # -2 * sum(log(c(0.9, 0.7, 0.6, 0.6, 0.8, 0.9))) = 3.6243815
  # Averaged over 6 sites = 0.6040636
  expect_equal(
    calc.deviance(siteData$obs, siteData$preds),
    0.6040636,
    tolerance = 1e-6
  )
  expect_equal(
    calc.deviance(siteData$obs, siteData$preds, calc.mean = FALSE),
    3.6243815,
    tolerance = 1e-6
  )
})

test_that("calc.deviance stops when obs and preds differ in length", {
  expect_error(
    calc.deviance(obs = c(1, 0, 1), preds = c(0.5, 0.5)),
    "equal length"
  )
})

# tweak.p ----------------------------------------------------------------------

test_that("tweak.p replaces 0 and 1 with the nearest in-between values", {
  expect_equal(tweak.p(c(0, 0.2, 0.8, 1)), c(0.2, 0.2, 0.8, 0.8))
})

test_that("tweak.p leaves values strictly between 0 and 1 unchanged", {
  probabilities <- c(0.1, 0.5, 0.9)

  expect_equal(tweak.p(probabilities), probabilities)
})

# calcSiteWeights --------------------------------------------------------------

test_that("calcSiteWeights down-weights background to balance presences", {
  response <- c(1, 1, 0, 0, 0, 0)

  weights <- calcSiteWeights(response)

  # 2 presences / 4 background sites = 0.5 for each background site
  expect_equal(weights, c(1, 1, 0.5, 0.5, 0.5, 0.5))
  expect_equal(sum(weights[response == 0]), sum(response == 1))
})

test_that("calcSiteWeights gives every site weight 1 when all are presences", {
  expect_equal(calcSiteWeights(c(1, 1, 1)), c(1, 1, 1))
})

test_that("calcSiteWeights returns NA for sites with a missing response", {
  response <- c(1, 0, 0, NA)

  # NA is left out of the counts: 1 presence / 2 background sites = 0.5
  expect_equal(calcSiteWeights(response), c(1, 0.5, 0.5, NA))
})

# boyceIndex -------------------------------------------------------------------

# With bins = 5, the bins are 0-0.2, 0.2-0.4, 0.4-0.6, 0.6-0.8 and 0.8-1.
# One prediction in each bin means every bin is equally available.
evenPredictions <- c(0.1, 0.3, 0.5, 0.7, 0.9)

test_that("boyceIndex is 1 when presences increase with predicted value", {
  # 0, 1, 2, 3 and 4 presences in the five bins
  presPredictions <- c(0.3, 0.5, 0.5, 0.7, 0.7, 0.7, 0.9, 0.9, 0.9, 0.9)

  expect_equal(boyceIndex(evenPredictions, presPredictions, bins = 5), 1)
})

test_that("boyceIndex is -1 when presences decrease with predicted value", {
  # 4, 3, 2, 1 and 0 presences in the five bins
  presPredictions <- c(0.1, 0.1, 0.1, 0.1, 0.3, 0.3, 0.3, 0.5, 0.5, 0.7)

  expect_equal(boyceIndex(evenPredictions, presPredictions, bins = 5), -1)
})

test_that("boyceIndex is NA when predictions fall in fewer than 3 bins", {
  allPredictions <- c(0.1, 0.1, 0.9, 0.9)

  expect_identical(
    boyceIndex(allPredictions, presPredictions = 0.9, bins = 5),
    NA_real_
  )
})

# safe_rbind -------------------------------------------------------------------

test_that("safe_rbind fills columns missing from the new row with NA", {
  existingRows <- data.frame(Name = "a", Value = 1)

  combined <- safe_rbind(existingRows, data.frame(Name = "b"))

  expect_equal(nrow(combined), 2)
  expect_equal(names(combined), c("Name", "Value"))
  expect_true(is.na(combined$Value[2]))
})

test_that("safe_rbind adds columns that only exist in the new row", {
  existingRows <- data.frame(Name = "a")

  combined <- safe_rbind(existingRows, list(Name = "b", Extra = 5))

  expect_equal(names(combined), c("Name", "Extra"))
  expect_equal(combined$Extra, c(NA, 5))
})
