## Unit tests for src/05-training-testing-data-prep-functions.R -- pure
## data.frame splitting logic that doesn't require a live SyncroSim scenario.

test_that("testTrainSplit splits each response level near the requested proportion", {
  source_wisdm_r("05-training-testing-data-prep-functions.R")

  set.seed(1)
  n <- 220
  dat <- data.frame(Response = rep(c(0, 1), each = n / 2), x1 = rnorm(n))

  result <- testTrainSplit(dat, trainProp = 0.75)

  expect_equal(nrow(result), n)
  expect_setequal(unique(result$UseInModelEvaluation), c(TRUE, FALSE))

  testFrac <- tapply(result$UseInModelEvaluation, result$Response, mean)
  expect_equal(as.numeric(testFrac), c(0.25, 0.25), tolerance = 0.05)
})

test_that("testTrainSplit errors when background (nodataValue) rows are present", {
  source_wisdm_r("05-training-testing-data-prep-functions.R")

  # NOTE: this documents a pre-existing bug, not desired behavior. `bg.dat` is
  # captured before `UseInModelEvaluation` is added to `dat`, so the final
  # `names(bg.dat) <- names(dat)` fails on a column-count mismatch whenever
  # nodataValue rows reach this function -- contradicting the function's own
  # doc comment ("Background points are ... read in and written out"). It's
  # unreachable in the real pipeline today because
  # 5-prepare-training-testing-data.R already strips nodataValue rows before
  # calling testTrainSplit(). If that's fixed upstream, update this test to
  # assert the background row is passed through untouched instead.
  set.seed(1)
  n <- 200
  dat <- data.frame(Response = rep(c(0, 1), each = n / 2), x1 = rnorm(n))
  dat <- rbind(dat, data.frame(Response = nodataValue, x1 = 0))

  expect_error(testTrainSplit(dat, trainProp = 0.8), "must be the same length")
})

test_that("testTrainSplit rejects an out-of-range training proportion", {
  source_wisdm_r("05-training-testing-data-prep-functions.R")
  dat <- data.frame(Response = rep(0:1, 100), x1 = rnorm(200))
  expect_error(testTrainSplit(dat, trainProp = 1.5), "Training Proportion")
})

test_that("testTrainSplit rejects fewer than 100 observations", {
  source_wisdm_r("05-training-testing-data-prep-functions.R")
  dat <- data.frame(Response = rep(0:1, 10), x1 = rnorm(20))
  expect_error(testTrainSplit(dat, trainProp = 0.7), "not advisable")
})

test_that("crossValidationSplit assigns every row to exactly one, roughly-balanced fold", {
  source_wisdm_r("05-training-testing-data-prep-functions.R")

  set.seed(1)
  n <- 500
  dat <- data.frame(
    Response = rep(c(0, 1), each = n / 2),
    UseInModelEvaluation = FALSE,
    ModelSelectionSplit = NA_integer_
  )

  result <- crossValidationSplit(dat, factorVars = NULL, nFolds = 5)

  expect_equal(nrow(result), n)
  expect_false(anyNA(result$ModelSelectionSplit))
  expect_setequal(unique(result$ModelSelectionSplit), 1:5)
  expect_true(all(abs(table(result$ModelSelectionSplit) - n / 5) <= 1))
})

test_that("crossValidationSplit rejects an invalid number of folds", {
  source_wisdm_r("05-training-testing-data-prep-functions.R")
  dat <- data.frame(
    Response = rep(0:1, 50),
    UseInModelEvaluation = FALSE,
    ModelSelectionSplit = NA_integer_
  )
  expect_error(crossValidationSplit(dat, factorVars = NULL, nFolds = 1), "greater than 1")
})
