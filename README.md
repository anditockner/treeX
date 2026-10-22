# treeX <img src="https://raw.githubusercontent.com/anditockner/treeX/refs/heads/main/man/figures/treeXLogo_500.png" align="right" style="width: 250px; min-width: 200px;"/>


**R-package for individual tree processing from ground-based point clouds.**

The functions in treeX were developed to analyze ground-based LiDAR data (PLS, MLS, TLS) from forest plots. 
It provides a modular framework for ground removal (including the generation of vertical slices), individual tree detection (via density clustering), instance segmentation (via iterative region-growing), and individual tree measurement (crown features are derived via alpha-hulls). 

Automatically derived parameters for diameter at breast height, tree height and crown basal height were benchmarked with conventional measurements in forest inventory (with a root mean square error of less than 2 cm dbh and 2 m height). See the citation section for more information.  

The package and all its functions were developed at the Institute of Forest Growth (WAFO), BOKU University, Vienna. 

A sample workflow can be found in the script "WORKFLOW_sample.R". 


# Installation

Before first use, make sure to also run the script "_install_for_treeX.R" once for installing the required libraries. 

       source("https://raw.githubusercontent.com/anditockner/treeX/refs/heads/main/_install_for_treeX.R")


After getting all dependencies, install the package from GitHub via

       pak::pak("anditockner/treeX")



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
