## -------------------------
## wisdm - R environment smoke tests
## ApexRMS, October 2026
## -------------------------

# Checks that the wisdm R conda environment matches src/wisdm-r-conda.yml
# (every pinned package loads, at the pinned version) and that the spatial
# libraries GDAL and PROJ work. Run these after changing a version in the yml
# file. They are skipped outside the wisdm R conda environment.

# Read the yml file ------------------------------------------------------------

condaLines <- readLines(file.path(srcDir, "wisdm-r-conda.yml"))

envName <- sub("^name:\\s*", "", grep("^name:", condaLines, value = TRUE))

# Lines like "  - r-terra=1.5_21" become condaName "r-terra", version "1.5_21"
pinLines <- sub("^\\s*-\\s*", "", grep("^\\s*-\\s*\\S+=", condaLines,
                                       value = TRUE))
pins <- data.frame(
  condaName = sub("=.*", "", pinLines),
  version = sub(".*=", "", pinLines)
)

# conda writes versions like "1.3_5"; R writes the same version as "1.3.5"
normalizeVersion <- function(version) {
  return(as.character(package_version(gsub("_", "-", version))))
}

skipUnlessInWisdmEnv <- function() {
  skip_if_not(
    basename(Sys.getenv("CONDA_PREFIX")) == envName,
    paste("only runs inside the", envName, "conda environment")
  )
}

# Tests ------------------------------------------------------------------------

test_that("R is the pinned version", {
  skipUnlessInWisdmEnv()

  pinnedVersion <- pins$version[pins$condaName == "r-base"]

  expect_identical(as.character(getRversion()), pinnedVersion)
})

test_that("every pinned R package loads at the pinned version", {
  skipUnlessInWisdmEnv()

  installedPackages <- .packages(all.available = TRUE)
  packagePins <- pins[pins$condaName != "r-base", ]

  for (i in seq_len(nrow(packagePins))) {
    condaName <- packagePins$condaName[i]

    # conda names are lowercase with an "r-" prefix; R names can have capitals
    # (e.g. r-prroc is the R package PRROC)
    packageName <- installedPackages[
      tolower(installedPackages) == sub("^r-", "", condaName)
    ][1]

    if (is.na(packageName)) {
      fail(paste(condaName, "is not installed"))
      next
    }

    expect_true(
      requireNamespace(packageName, quietly = TRUE),
      label = paste(packageName, "loads")
    )
    expect_identical(
      as.character(utils::packageVersion(packageName)),
      normalizeVersion(packagePins$version[i]),
      label = paste(packageName, "version")
    )
  }
})

test_that("terra can write and read a GeoTIFF (GDAL works)", {
  skipUnlessInWisdmEnv()

  rasterFile <- tempfile(fileext = ".tif")
  original <- terra::rast(
    nrows = 10, ncols = 10,
    xmin = 0, xmax = 10, ymin = 0, ymax = 10,
    crs = "EPSG:4326",
    vals = 1:100
  )

  terra::writeRaster(original, rasterFile)
  reloaded <- terra::rast(rasterFile)

  expect_equal(dim(reloaded), c(10, 10, 1))
  expect_equal(terra::global(reloaded, "sum")[[1]], sum(1:100))
  expect_true(nzchar(terra::gdal()))
})

test_that("sf can reproject a point (PROJ works)", {
  skipUnlessInWisdmEnv()

  point <- sf::st_sfc(sf::st_point(c(-120, 40)), crs = 4326)

  projected <- sf::st_transform(point, 3857)

  # Web Mercator x = longitude * 20037508.34 / 180 = -13358338.9
  expect_equal(
    unname(sf::st_coordinates(projected)[1, "X"]),
    -13358338.9,
    tolerance = 1e-6
  )
})
