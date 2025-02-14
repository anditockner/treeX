R-package for tree detection and instance segmentation from ground-based point clouds. 

A sample workflow can be found in the script "WORKFLOW_sample.R". 

Installing the package from GitHub: 
devtools::install_github("https://github.com/anditockner/treeX")

The software generates a folder structure given in the path specified in dirPath. 
All files will be renamed according to a file identifier specified in fileFinder. 
The input file is specified in LASfile. 


Main functions: 

# preparing folder structure with pre-processing of point cloud and creating ground model
extractVegetation(LASfile, fileFinder)

# detecting tree positions and measuring diameter at 1.30 m (DBH)
clustSplit(fileFinder = fileFinder)

# instance segmentation starting from list defining tree DBH positions
crownFeel(fileFinder)

# calculating individual tree metrics (height, crown projection area, crown hull volume)
computeTreeParams(fileFinder)
