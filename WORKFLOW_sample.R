#devtools::install_github("anditockner/treeX")  
library(treeX)

# files and folder structure will be saved here
# (use a path that exists on your system, e.g. "D:/try/" on Windows
#  or "~/treeX_try/" on macOS/Linux)
dirPath <- path.expand("~/treeX_try/")
if(!dir.exists(dirPath)) dir.create(dirPath)
setwd(dirPath)

#specify the input file
LASfile <- list.files(dirPath, pattern = "laz", full.name = T)

# a characteristic name to save all output files
fileFinder <- "spzaiil"

voxelSize <- 20
# edge length of voxel for crown segmentation
# -> 20 cm to speed up the process
#    usually far too rough, decrease it down to 5 (2) cm for better results

# ground detection
try(extractVegetation(LASfile, fileFinder))

# individual tree detection
try(clustSplit(fileFinder = fileFinder))

# instance segmentation, region growing
try(crownFeel(fileFinder, voxelSize = voxelSize, tileClipping = 0,
              retainPointClouds = T))

# measure parameters (height, crown basal height, hull volume)
try(computeTreeParams(fileFinder,
                      voxelSize = voxelSize,
                      crownParameters = T))


