## Entrypoint for the wisdm R test suite.
## Run from the repo root with: Rscript tests/R/testthat.R

library(testthat)

find_wisdm_repo_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "src", "package.xml"))) {
      return(path)
    }
    parent <- dirname(path)
    if (parent == path) {
      stop("Could not locate the wisdm repo root (src/package.xml) above ", start)
    }
    path <- parent
  }
}

repoRoot <- find_wisdm_repo_root()
test_dir(file.path(repoRoot, "tests", "R", "testthat"))
