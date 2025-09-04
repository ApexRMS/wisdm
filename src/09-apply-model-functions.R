## ------------------------------
## wisdm - apply model functions  
## ApexRMS, April 2022       
## ------------------------------ 

# Update breakpoint function ---------------------------------------------------
# Function to time code by returning a clean string of time since this function was last called

updateBreakpoint <- function() { 

  # Calculate time since last breakpoint
  newBreakPoint <- proc.time()
  elapsed <- (newBreakPoint - currentBreakPoint)['elapsed']
  
  # Update current breakpoint
  currentBreakPoint <<- newBreakPoint
  
  # Return cleaned elapsed time
  if (elapsed < 60) {
    return(paste0(round(elapsed), "sec"))
  } else if (elapsed < 60^2) {
    return(paste0(round(elapsed / 60, 1), "min"))
  } else
    return(paste0(round(elapsed / 60 / 60, 1), "hr"))
}

# detect system memory function ------------------------------------------------
# Cross-platform system memory detection (GB)

detect_system_memory <- function() {
  out <- list(total_gb = NA_real_, available_gb = NA_real_, method = "unknown")
  # Preferred: ps package (works on Windows/macOS/Linux)
  if (requireNamespace("ps", quietly = TRUE)) {
    m <- ps::ps_system_memory()
    out$total_gb     <- as.numeric(m$total)     / (1024^3)
    out$available_gb <- as.numeric(m$available) / (1024^3)
    out$method <- "ps::ps_system_memory()"
    return(out)
  }
  sys <- tolower(Sys.info()[["sysname"]])

  if (sys == "linux") {
    # Read /proc/meminfo (kB)
    mi <- readLines("/proc/meminfo", warn = FALSE)
    num <- function(key) {
      line <- mi[grepl(paste0("^", key, ":"), mi)]
      if (!length(line)) return(NA_real_)
      as.numeric(sub(".*?([0-9]+).*", "\\1", line)) / (1024^2)  # -> GB
    }
    out$total_gb     <- num("MemTotal")
    out$available_gb <- num("MemAvailable")
    out$method <- "/proc/meminfo"
    return(out)
  }

  if (sys == "windows") {
    # Try PowerShell CIM (bytes)
    ps_cmd <- 'powershell -NoProfile -Command "(Get-CimInstance Win32_OperatingSystem | Select-Object TotalVisibleMemorySize,FreePhysicalMemory | ConvertTo-Json -Compress)"'
    j <- tryCatch(system(ps_cmd, intern = TRUE), error = function(e) "")
    if (nzchar(paste(j, collapse = ""))) {
      dat <- tryCatch(jsonlite::fromJSON(paste(j, collapse = "")), error = function(e) NULL)
      if (!is.null(dat$TotalVisibleMemorySize)) {
        out$total_gb     <- as.numeric(dat$TotalVisibleMemorySize) / (1024^2) # kB -> GB
        out$available_gb <- as.numeric(dat$FreePhysicalMemory)     / (1024^2)
        out$method <- "PowerShell CIM"
        return(out)
      }
    }
    # Fallback: wmic (kB) — deprecated but often present
    w <- tryCatch(system('wmic OS get TotalVisibleMemorySize,FreePhysicalMemory /Value', intern = TRUE), error = function(e) "")
    if (length(w)) {
      getv <- function(key) {
        line <- w[grepl(paste0("^", key, "="), w)]
        if (!length(line)) return(NA_real_)
        as.numeric(sub(".*=(\\d+).*", "\\1", line)) / (1024^2)  # kB -> GB
      }
      out$total_gb     <- getv("TotalVisibleMemorySize")
      out$available_gb <- getv("FreePhysicalMemory")
      out$method <- "wmic"
      return(out)
    }
  }

  if (sys == "darwin") {
    # macOS: total via sysctl (bytes)
    total_b <- suppressWarnings(as.numeric(system("sysctl -n hw.memsize", intern = TRUE)))
    # Available (approx): free + inactive pages * page size
    pagesize <- suppressWarnings(as.numeric(system("sysctl -n hw.pagesize", intern = TRUE)))
    vm <- tryCatch(system("vm_stat", intern = TRUE), error = function(e) character())
    getp <- function(pattern) {
      line <- vm[grepl(pattern, vm)]
      if (!length(line)) return(NA_real_)
      as.numeric(sub(".*:\\s*([0-9]+)\\.", "\\1", line))
    }
    free_p     <- getp("^Pages free")
    inactive_p <- getp("^Pages inactive")
    avail_b <- if (is.finite(pagesize) && is.finite(free_p) && is.finite(inactive_p))
      (free_p + inactive_p) * pagesize else NA_real_
    out$total_gb     <- if (is.finite(total_b)) total_b / (1024^3) else NA_real_
    out$available_gb <- if (is.finite(avail_b)) avail_b / (1024^3) else NA_real_
    out$method <- "sysctl/vm_stat"
    return(out)
  }

  out
}

# Session setup function -------------------------------------------------------
# Sets up terra, GDAL, BLAS/MKL env vars

setup_session <- function(
  ssim_temp_dir,                 # where to put temps (your TransferDirectory)
  concurrent_sessions = 1L,      # how many R sessions on this host
  memmax_gb          = NULL,     # hard cap per session (GB). If NULL, inferred.
  total_ram_gb       = NULL,     # optional: total machine RAM (GB) to infer memmax
  gdal_cache_mb      = NULL,     # per-process GDAL cache (MB). If NULL, inferred.
  block_rows         = NULL,     # residuals rows per block. If NULL, inferred.
  desync_max_sec     = 0,        # random delay [0, desync_max_sec] to stagger I/O
  threads_per_lib    = 1L        # OMP/BLAS/MKL threads per process
) {

  stopifnot(dir.exists(ssim_temp_dir))
  S <- as.integer(concurrent_sessions)

  # Choose memmax (GB) ---
  if (is.null(memmax_gb)) {
    if (is.numeric(total_ram_gb) && is.finite(total_ram_gb) && total_ram_gb > 0) {
      # keep aggregate usage ≲ ~60% of RAM across sessions
      memmax_gb <- max(0.5, min(8, 0.6 * total_ram_gb / S))
    } else {
      memmax_gb <- if (S >= 16) 1.0 else if (S >= 8) 1.5 else if (S >= 4) 2.0 else 3.0
    }
  }

  # GDAL cache (MB): ~ memmax/4, clamped [128, 1024] unless provided ---
  if (is.null(gdal_cache_mb)) {
    gdal_cache_mb <- max(128L, min(1024L, as.integer(round(memmax_gb * 256))))
  }

  # Residuals block size ---
  if (is.null(block_rows)) {
    block_rows <- if (memmax_gb <= 1.0) 256L else if (memmax_gb <= 2.0) 384L else 512L
  }
  block_rows <- max(64L, as.integer(block_rows))

  # Apply to terra & GDAL/BLAS ---
  terraOptions(tempdir = ssim_temp_dir, memmax = memmax_gb)
  Sys.setenv(
    GDAL_CACHEMAX        = as.character(gdal_cache_mb),
    OMP_NUM_THREADS      = as.character(threads_per_lib),
    OPENBLAS_NUM_THREADS = as.character(threads_per_lib),
    MKL_NUM_THREADS      = as.character(threads_per_lib)
  )

  # --- Optional: stagger starts to reduce I/O spikes ---
  if (desync_max_sec > 0) Sys.sleep(stats::runif(1, 0, desync_max_sec))

  # Return knobs the rest of the script uses
  list(
    MEMMAX_GB          = memmax_gb,
    GDAL_CACHE_MB      = gdal_cache_mb,
    BLOCK_ROWS         = block_rows
  )
}

# Build cov stack function -----------------------------------------------------
# Helper: build & align a covariate stack for requested variables (v1.5-21 safe)

build_cov_stack <- function(cov_data, var_names, factor_vars) {
  
  rows  <- cov_data[cov_data$CovariatesID %in% var_names, , drop = FALSE]
  rows  <- rows[!duplicated(rows$CovariatesID), , drop = FALSE]
  paths <- rows$RasterFilePath
  
  if (any(!file.exists(paths))) {
    missingTifs <- rows$CovariatesID[!file.exists(paths)]
    stop("Missing rasters in Covariate Data: ", paste(missingTifs, collapse = ", "))
  }
  if (sum(file.access(paths), mode = 0) != 0) {
    stop("Unreadable raster: ", paths[(file.access(paths) != 0)][1])
  }
  
  s <- rast(paths)
  names(s) <- rows$CovariatesID
  
  for (nm in intersect(factor_vars, names(s))){
    s[[nm]] <- as.factor(s[[nm]])
  } 
  return(s)
}

# --------------------------------------------------------------------
# Probability prediction block function (single-pass INT output)
# - Coerces factors to precomputed levels
# - Applies restriction in the same pass (if present)
# - Returns integer 0..100 for direct INT2S writing
# --------------------------------------------------------------------

predict_block_int <- function(model, data, factor_levels, predictFct, restrict = NULL, ...) {

  # 1) Force base data.frame with atomic vectors (no S4/list columns)
  if (!is.data.frame(data)){ data <- as.data.frame(data, stringsAsFactors = FALSE)}
  for (nm in names(data)) {
    if (methods::is(data[[nm]], "S4")) {
      data[[nm]] <- as.vector(data[[nm]])
    }
  }

  # 2) Apply factor levels (only if the column exists)
  if (length(factor_levels)) {
    for (j in names(factor_levels)) if (j %in% names(data)) {
      lv <- factor_levels[[j]]
      data[[j]] <- if (is.null(lv)) factor(data[[j]]) else factor(data[[j]], levels = lv)
    }
  }
  
  # 3) Apply prediction function
  p <- predictFct(model = model, x = data)  # numeric [0,1]
  
  # 4) Apply restriction (if present)
  if (!is.null(restrict)) {
    restr <- as.vector(restrict)
    p[is.na(restr)] <- NA
    if (is.numeric(restr)) { p[restr == 0] <- NA }
  }

  # 5) Convert to integer [0,100] (probability as percent)
  as.integer(round(p * 100))
}

# --------------------------------------------------------------------
# Fast MESS/MoD helpers (block friendly)
# --------------------------------------------------------------------
prepare_mess <- function(train.dat, var_names) {
  td <- as.data.frame(train.dat[, var_names, drop = FALSE])
  for (nm in names(td)) if (is.factor(td[[nm]])) td[[nm]] <- as.numeric(td[[nm]])
  mins <- vapply(td, min, numeric(1), na.rm = TRUE)
  maxs <- vapply(td, max, numeric(1), na.rm = TRUE)
  n    <- nrow(td)
  sorted_cols <- lapply(td, function(col) sort(col))
  list(var_names = var_names, mins = mins, maxs = maxs, n = n, sorted_cols = sorted_cols)
}

range_proximity_pct <- function(x_mat, mins, maxs) {
  span  <- matrix(rep(maxs - mins, each = nrow(x_mat)), nrow = nrow(x_mat))
  span[span == 0] <- NA_real_
  left  <- sweep(x_mat, 2, mins, `-`)
  right <- sweep(matrix(rep(maxs, each = nrow(x_mat)), nrow = nrow(x_mat)), 2, x_mat, `-`)
  pmin(left, right) / span * 100
}

percent_rank_distance <- function(x_mat, sorted_cols, n) {
  mats <- lapply(seq_along(sorted_cols), function(j) {
    col_sorted <- sorted_cols[[j]]
    x          <- x_mat[, j]
    countLess  <- findInterval(x, col_sorted, left.open = FALSE, rightmost.closed = TRUE)
    p <- (countLess / n) * 100
    2 * pmin(p, 100 - p)
  })
  do.call(cbind, mats)
}

CalcMESS_fast <- function(rast_block, prep) {
  x <- as.data.frame(rast_block)[, prep$var_names, drop = FALSE]
  for (nm in names(x)) if (is.factor(x[[nm]])) x[[nm]] <- as.numeric(x[[nm]])
  x_mat <- as.matrix(x)
  row_has_na <- !stats::complete.cases(x_mat)
  
  prox_mat     <- range_proximity_pct(x_mat, prep$mins, prep$maxs)
  prox_mat_neg <- -prox_mat
  min_idx      <- max.col(prox_mat_neg, ties.method = "first")
  min_val      <- prox_mat[cbind(seq_len(nrow(x_mat)), min_idx)]
  
  inside <- !row_has_na & is.finite(min_val) & (min_val > 0)
  if (any(inside)) {
    rank_mat      <- percent_rank_distance(x_mat[inside, , drop = FALSE], prep$sorted_cols, prep$n)
    rank_mat_neg  <- -rank_mat
    min_idx_inside<- max.col(rank_mat_neg, ties.method = "first")
    min_val_inside<- rank_mat[cbind(seq_len(nrow(rank_mat)), min_idx_inside)]
    min_idx[inside] <- min_idx_inside
    min_val[inside] <- min_val_inside
  }
  
  idx_names <- prep$var_names[min_idx]
  if (any(row_has_na)) {
    idx_names[row_has_na] <- NA_character_
    min_val[row_has_na]   <- NA_real_
  }
  data.frame(Indx = idx_names, output = min_val, row.names = NULL)
}
