setwd("D:/treeX_git/")
dyn.load("edc.dll") #C-library for circle fits
source("all_initialize.R")
source("classify_ground.R")
source("detect_trees.R")
source("instance_segmentation.R")
source("parameter_measuring.R")


dirPath <- "D:/try/" # files and folder structure will be saved here
if(!dir.exists(dirPath)) dir.create(dirPath)
setwd(dirPath)


#specify the input file
LASfile <- list.files("D:/PointCloud_Processing/", pattern = "laz", full.name = T)[2]

# a characteristic name to save all output files
fileFinder <- "spzaiil"

voxelSize <- 20 
# edge length of voxel for crown segmentation -> 20 cm to speed up the process

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


