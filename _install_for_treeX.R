# ---------------------------------------------------------------------------
# Install script for treeX dependencies
# Supports Windows, macOS (Intel and Apple Silicon) and Linux.
# Run this once before installing treeX itself.
# ---------------------------------------------------------------------------

# macOS notes ---------------------------------------------------------------
# * treeX now compiles a small C++ routine from source; install the Apple
#   command line tools once if you have not already:   xcode-select --install
# * For 3D plots the 'rgl' package uses XQuartz: https://www.xquartz.org/
# * All CRAN packages below are available as macOS binaries, so no further
#   system libraries are required.

if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

is_windows <- .Platform$OS.type == "windows"
is_macos   <- Sys.info()[["sysname"]] == "Darwin"

# --- core dependencies (Depends/Imports of treeX) --------------------------
core <- c("alphahull", "alphashape3d", "conicfit", "dbscan", "dismo",
          "doParallel", "dplyr", "foreach", "lidR", "ks", "Metrics",
          "Morpho", "plotrix", "plyr", "RANN", "sf", "sp", "spatstat",
          "stringi", "VoxR", "data.table")

# --- optional helpers used by some functions (Suggests) --------------------
extra <- c("densityClust", "DescTools", "digest", "geosphere",
           "getPass", "igraph", "jpeg", "MASS", "matrixStats", "mgcv",
           "PBSmodelling", "png", "pracma", "raster", "rgl", "Rvcg",
           "RCSF", "flexclust", "benchmarkme")

# progress bar for parallel runs on macOS/Linux (not needed on Windows)
if (!is_windows) extra <- c(extra, "mcprogress")

install.packages(setdiff(c(core, extra), rownames(installed.packages())),
                 dependencies = TRUE)

# lidR from r-universe gets you the latest version (optional; CRAN also works)
# install.packages("lidR", repos = c("https://r-lidar.r-universe.dev",
#                                    "https://cloud.r-project.org"))

# --- optional: ADPclust (archived on CRAN, GitHub only) --------------------
# Loaded defensively by treeX when present, but none of its functions are
# actually called - skipping it is fine.
if (!requireNamespace("ADPclust", quietly = TRUE)) {
  tryCatch(pak::pak("ethanyxu/ADPclust"),
           error = function(e) {
             message("ADPclust could not be installed - that is fine, ",
                     "treeX does not require it.")
           })
}

# --- optional: TreeLS (archived on CRAN, GitHub only) ----------------------
# treeX no longer requires TreeLS: voxel thinning falls back to lidR when
# TreeLS is absent. Installing it is optional. On macOS its compilation
# needs the CRAN gfortran toolchain (https://mac.r-project.org/tools/) -
# without it, linking fails ("ld: library 'emutls_w' not found"); simply
# skip it in that case.
tryCatch(pak::pak("tiagodc/TreeLS"),
         error = function(e) {
           message("TreeLS could not be installed - that is fine, ",
                   "treeX will use its lidR-based fallback instead.")
         })

# --- finally: install treeX itself (compiles its C++ code from source) -----
# pak::pak("anditockner/treeX")
