R-package for tree detection and instance segmentation from ground-based point clouds. 

A sample workflow can be found in the script "WORKFLOW_sample.R". 

# Installation

Install the package from GitHub via

       devtools::install_github("https://github.com/anditockner/treeX")

Before first use, make sure to also run the script "_install_for_treeX.R" once for installing the required libraries. 

       source("https://raw.githubusercontent.com/anditockner/treeX/refs/heads/main/_install_for_treeX.R")


The software generates a folder structure given in the path specified in dirPath. 
All files will be renamed according to a file identifier specified in fileFinder. 
The input file is specified in LASfile. 


# Main Functions

* preparing folder structure with pre-processing of point cloud and creating ground model

       extractVegetation(LASfile, fileFinder)

* detecting tree positions and measuring diameter at 1.30 m (DBH)

       clustSplit(fileFinder)

* instance segmentation starting from list defining tree DBH positions

       crownFeel(fileFinder)

* calculating individual tree metrics (height, crown projection area, crown hull volume)

       computeTreeParams(fileFinder)



# Citation

For citing this methodology and for more details on the implementation, please refer to following paper

Tockner et al. 2022: Automatic tree crown segmentation using dense forest point clouds from Personal Laser Scanning (PLS)

https://doi.org/10.1016/j.jag.2022.103025
