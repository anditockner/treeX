# ---------------------------------------------------------------------------
# Cross-platform compatibility helpers (added for macOS support)
# ---------------------------------------------------------------------------

#' Register a parallel backend that works on every operating system
#'
#' On Windows a PSOCK cluster is created and registered via doParallel.
#' On macOS and Linux doParallel is registered in fork mode
#' (equivalent to the previously used doMC backend, but without the extra
#' dependency - doMC is not needed any more).
#'
#' @param nr_cores number of worker processes
#' @param outfile passed to parallel::makeCluster on Windows ("" shows worker
#'   output in the console)
#' @return the cluster object on Windows, otherwise NULL (fork backend)
#' @keywords internal
registerParallelBackend <- function(nr_cores, outfile = NULL) {
  if (.Platform$OS.type == "windows") {
    cl <- if (is.null(outfile)) {
      parallel::makeCluster(nr_cores)
    } else {
      parallel::makeCluster(nr_cores, outfile = outfile)
    }
    doParallel::registerDoParallel(cl)
    cat("   (on Windows using doParallel PSOCK cluster)\n")
    return(invisible(cl))
  } else {
    doParallel::registerDoParallel(cores = nr_cores)
    cat("   (on Unix-like system (macOS/Linux) using doParallel fork backend)\n")
    return(invisible(NULL))
  }
}

#' Attach packages only if they are installed
#'
#' Used for the defensive package loads inside functions that run on parallel
#' workers. Several of the historically loaded packages (e.g. recexcavAAR,
#' TreeLS) are no longer on CRAN and may be missing, especially on macOS;
#' a missing optional package must not abort the whole run.
#'
#' @param pkgs character vector of package names
#' @keywords internal
attachIfAvailable <- function(pkgs) {
  for (p in pkgs) {
    if (suppressWarnings(requireNamespace(p, quietly = TRUE))) {
      suppressPackageStartupMessages(
        library(p, character.only = TRUE)
      )
    }
  }
  invisible(NULL)
}

#' Voxel-based point cloud thinning (one random point per voxel)
#'
#' Drop-in replacement for TreeLS::tlsSample(las, TreeLS::smp.voxelize(res)).
#' TreeLS is archived and only installable from GitHub, which frequently
#' fails to compile on macOS; lidR provides the same thinning natively.
#' If TreeLS is installed it is still used, so results stay identical to
#' earlier versions of treeX.
#'
#' @param las a LAS object
#' @param res voxel edge length in point cloud units (usually meters)
#' @return the thinned LAS object
#' @keywords internal
voxelSample <- function(las, res) {
  if (requireNamespace("TreeLS", quietly = TRUE)) {
    return(TreeLS::tlsSample(las, TreeLS::smp.voxelize(res)))
  }
  lidR::decimate_points(las, lidR::random_per_voxel(res = res, n = 1L))
}

#' Locate the CloudCompare executable on the current platform
#'
#' @param userPath optional user-supplied path to the CloudCompare
#'   executable or its directory; if it exists it wins.
#' @return quoted command string for the CloudCompare binary
#' @keywords internal
findCloudCompare <- function(userPath = NULL) {
  # user override (either full binary path or its directory)
  if (!is.null(userPath) && nzchar(userPath)) {
    clean <- gsub('"', '', userPath)
    if (file.exists(clean) && !dir.exists(clean)) return(shQuote(clean))
    exe <- if (.Platform$OS.type == "windows") "CloudCompare.exe" else "CloudCompare"
    cand <- file.path(clean, exe)
    if (file.exists(cand)) return(shQuote(cand))
  }

  sysname <- Sys.info()[["sysname"]]
  candidates <- if (.Platform$OS.type == "windows") {
    c("C:/Program Files/CloudCompare/CloudCompare.exe",
      "C:/Program Files (x86)/CloudCompare/CloudCompare.exe")
  } else if (identical(sysname, "Darwin")) {
    c("/Applications/CloudCompare.app/Contents/MacOS/CloudCompare",
      file.path(Sys.getenv("HOME"),
                "Applications/CloudCompare.app/Contents/MacOS/CloudCompare"))
  } else {
    c(Sys.which("CloudCompare"), Sys.which("cloudcompare.CloudCompare"))
  }
  candidates <- candidates[nzchar(candidates)]
  hit <- candidates[file.exists(candidates)][1]
  if (!is.na(hit) && length(hit) == 1) return(shQuote(hit))

  stop("CloudCompare executable not found.\n",
       if (identical(sysname, "Darwin")) {
         "  On macOS install it from https://www.cloudcompare.org/ (or 'brew install --cask cloudcompare')\n  so that it is available at /Applications/CloudCompare.app.\n"
       } else if (.Platform$OS.type == "windows") {
         "  Expected under C:/Program Files/CloudCompare/.\n"
       } else {
         "  Make sure 'CloudCompare' is on the PATH (e.g. snap/flatpak wrapper).\n"
       },
       "  Alternatively pass the path explicitly (argument 'pathCloudCompare' / variable 'Path.CloudCompare').",
       call. = FALSE)
}

#' Run a system command, portably
#'
#' base::system(..., show.output.on.console=) is Windows-only; this wrapper
#' maps the intent (silent or not) onto the arguments each OS supports.
#'
#' @param command the command line to run
#' @param silent suppress command output
#' @param wait wait for the command to finish
#' @keywords internal
runCommand <- function(command, silent = FALSE, wait = TRUE) {
  if (.Platform$OS.type == "windows") {
    if (silent) {
      system(command, show.output.on.console = FALSE, wait = wait)
    } else {
      system(command, wait = wait)
    }
  } else {
    if (silent) {
      system(paste(command, ">/dev/null 2>&1"), wait = wait)
    } else {
      system(command, wait = wait)
    }
  }
}
