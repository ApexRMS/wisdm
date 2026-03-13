## -------------------------
## wisdm - constants
## ApexRMS, October 2025
## -------------------------

# Raster nodata sentinel value — must match nodataValue in setup_functions.py
nodataValue <- -9999

# Response value used to mark generated background/pseudo-absence points
backgroundValue <- -9998

# GAM smoothing methods
gamSmoothingMethodCW <- data.frame(codeTerm = c("ts", "cc", "re","gp"), displayTerm = c("Thin plate regression splines with shrinkage-to-zero",  "Cyclic cubic regression splines", "Random effects", "Gaussian process smooths"))
