
# this script is to be run on cmd and uses 4 arguments:
#
# [1] path to .las or .laz input ground-based point cloud
# [2] path to configuration file (config.R) defining hyper-parameters
# [3] (optional) path for saving output files
# [4] (optional) string for renaming output files (fileFinder)

# sample line to run it on cmd: 
# (use R 4.3.2 or 4.3.3)
# "C:/Programme/R/R-4.5.2/bin/x64/Rscript" 
#        "D:/WORKFLOW_treeX_cmd.R" "D:/input.laz" "D:/config_benchmark.R" "D:/output/" "firstSet"
# path to R + script workflow           [1]               [2]                [3]          [4]

time1 <- Sys.time()
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  print("Please specify .las or .laz file to be processed!")
  stop("Requires command line argument.")
}

# [1] input point cloud
LASfile <- args[1]
if (!file.exists(LASfile)) {
  print(paste0("File ", LASfile, " not found... Terminating."))
  stop("Invalid LASfile specified.")
}

# [2] config.R
if (length(args) >= 2) {
  if (!file.exists(args[2])) {
    print(paste0("Config file ", args[2], " not found... Terminating."))
    stop("Invalid config-file specified.")
  }
  source(args[2])
}

# [3] directory
dirPath <- dirname(LASfile) # files and folder structure will be created
if (length(args) >= 3) {
  dirPath <- args[3]
}

# [3] name
fileFinder <- gsub(".laz", "", basename(LASfile))
fileFinder <- gsub(".las", "", fileFinder)
if (length(args) >= 4) {
  fileFinder <- args[4]
}

cat("Loading treeX package... ")
so <- tryCatch(find.package("treeX"),
  error = function(error_condition) {
    cat("(installing from github)")
    devtools::install_github(
      "https://github.com/anditockner/treeX",
      upgrade = "never"
    )
  }
)
library("treeX")
cat("done!\n")

if (!dir.exists(dirPath)) dir.create(dirPath)
setwd(dirPath) # all the files and folder structure will be created here
cat("Working directory is:", dirPath, "\n")

cat("Starting treeX workflow\n")

cat(paste0("LASfile = ", LASfile), "\n")
cat(paste0("fileFinder = ", fileFinder), "\n")
cat("\n\n\n\n")

try(extractVegetation(LASfile, fileFinder,
  selector = selector,
  clip.radius = clip.radius, 
  clip.trajectory.distance = clip.trajectory.distance,
  exportSlice.upperLimit = exportSlice.upperLimit,
  exportSlice.lowerLimit = exportSlice.lowerLimit
))

try(clustSplit(
  fileFinder = fileFinder, 
  filterINT = filterINT,
  nr_cores = nr_cores_tree_detection,
  retainPointClouds = T
))

try(crownFeel(fileFinder,
  voxelSize = voxelSize, 
  limitStems = limitStems,
  limitShare = limitShare, 
  zScale = zScale,
  selector = selector,
  tileClipping = 2, 
  retainPointClouds = F
))

#
# try(computeTreeParams(fileFinder, voxelSize = voxelSize,
#                      limitShare = limitShare, zScale = zScale, nr_cores = 2,
#                      writeLAZ = F, retainPointClouds = F, crownParameters = F))
#

{
  time2 <- Sys.time()
  timeDiff <- round(time2 - time1, 2)
  cat("\n\n\nEnding script on", format(Sys.time()), "\n")
  cat("All jobs done in a ")
  print.difftime(timeDiff)
  cat("\n\n\n")
}
