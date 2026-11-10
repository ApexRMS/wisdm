# How testthat works

[testthat](https://testthat.r-lib.org/) is the most widely used testing
package for R. You write small blocks of code that call a function and check
the result. testthat runs them all and reports which passed and which failed.

For when and where to add tests in wisdm, see [../README.md](../README.md).

## What happens when you run the tests

`Rscript tests/r/run-tests.R` calls `testthat::test_dir("tests/r")`, which:

1. Sets the working directory to `tests/r`.
2. Runs every file whose name starts with `setup` (here,
   `setup-source-helpers.R` loads the wisdm functions from `src/`).
3. Runs every file whose name starts with `test-`, one `test_that()` block at
   a time.
4. Prints a summary.

Files that don't start with `setup`, `helper` or `test` (like `run-tests.R` or
this file) are ignored.

## A test, piece by piece

```r
test_that("tweak.p replaces 0 and 1 with the nearest in-between values", {
  probabilities <- c(0, 0.2, 0.8, 1)

  tweaked <- tweak.p(probabilities)

  expect_equal(tweaked, c(0.2, 0.2, 0.8, 0.8))
})
```

- **`test_that("description", { ... })`** groups the code for one behaviour.
  The description is shown when the test fails, so write it as a sentence
  that says what *should* happen.
- **Set up, run, check.** Most tests follow this pattern: build a small
  input, call the function, then check the result.
- **`expect_*()`** functions do the checking. If the check fails, testthat
  records a failure and moves on to the next check.

A `test_that()` block can contain more than one `expect_*()`, but keep each
block focused on one behaviour.

## Common expectations

| Expectation | Passes when |
|---|---|
| `expect_equal(x, y)` | `x` and `y` are equal (allows tiny rounding differences) |
| `expect_identical(x, y)` | `x` and `y` are exactly the same, including type |
| `expect_true(x)` / `expect_false(x)` | `x` is `TRUE` / `FALSE` |
| `expect_length(x, n)` | `x` has `n` elements |
| `expect_type(x, "double")` | `x` has that basic type |
| `expect_s3_class(x, "data.frame")` | `x` is (or inherits from) that class |
| `expect_named(x, c("a", "b"))` | `x` has those names |
| `expect_gt(x, y)` / `expect_lt(x, y)` | `x` is greater / less than `y` |
| `expect_match(x, "pattern")` | the string `x` matches the pattern |
| `expect_error(code, "message")` | `code` throws an error containing the message |
| `expect_warning(code, "message")` | `code` gives a warning containing the message |
| `expect_no_error(code)` | `code` runs without an error |

The full list is in the
[testthat reference](https://testthat.r-lib.org/reference/index.html).

Use `expect_equal()` for decimals, since computers store them with tiny
rounding errors. If your expected value is rounded (e.g. to 7 decimal
places), add a tolerance: `expect_equal(x, 0.6040636, tolerance = 1e-6)`.

### Checking for errors and warnings

Wrap the code that should fail in the expectation. The second argument is
part of the expected message:

```r
test_that("calc.deviance stops when obs and preds differ in length", {
  expect_error(
    calc.deviance(obs = c(1, 0, 1), preds = c(0.5, 0.5)),
    "equal length"
  )
})
```

### Using a fixture file

Use `test_path()` to find files in `tests/r`, whatever folder you run from:

```r
siteData <- read.csv(test_path("fixtures", "deviance-example.csv"))
```

If a test needs a package that might not be installed, start it with
`skip_if_not_installed("mgcv")` so it is skipped instead of failing.

## Running tests

From the repository root:

| To run | Command |
|---|---|
| All R tests | `Rscript tests/r/run-tests.R` |
| One file (from an R console) | `testthat::test_file("tests/r/test-00-helper-functions.R")` |

The summary line shows the number of failures, warnings, skips and passes,
e.g. `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 10 ]`. When a check fails, testthat
shows the file, the line number, the test
description, and how the values differ:

```
── Failure (test-00-helper-functions.R:41:3): calcSiteWeights down-weights ... ──
`weights` not equal to c(1, 1, 0.25, 0.25, 0.25, 0.25).
4/6 mismatches (average diff: 0.25)
```

## Going further

Features not covered here, for when you need them. wisdm's R environment is
pinned to R 4.1, so it may have an older testthat. Check the version notes
on each page.

- [Reference: all functions](https://testthat.r-lib.org/reference/index.html)
- [Special files](https://testthat.r-lib.org/articles/special-files.html):
  setup, helper and teardown files
- [Test fixtures](https://testthat.r-lib.org/articles/test-fixtures.html):
  temporary files, options and clean-up with `withr`
- [Mocking](https://testthat.r-lib.org/articles/mocking.html): replace a
  function (e.g. a SyncroSim call) with a fake one during a test
- [Snapshot tests](https://testthat.r-lib.org/articles/snapshotting.html):
  compare printed output, plots or files against a saved copy
- [Skipping tests](https://testthat.r-lib.org/articles/skipping.html)
- [Testing challenging functions](https://testthat.r-lib.org/articles/challenging-tests.html):
  randomness, graphics, user interaction
- [testthat 3rd edition](https://testthat.r-lib.org/articles/third-edition.html):
  what changed and how to opt in
- [R Packages, "Testing basics" chapter](https://r-pkgs.org/testing-basics.html):
  a longer, beginner-friendly introduction (written for R packages, so its
  folder names and `devtools` commands differ from wisdm's setup)
