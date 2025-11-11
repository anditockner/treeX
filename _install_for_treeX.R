install.packages("devtools", repos = "https://cran.rstudio.com/", dep = TRUE)

library("devtools")
{
  as <- find_rtools()
  if(!as){
    stop("\nRtools not found! It is required for installing TreeLS!\n\n")
  }
}

devtools::install_github('tiagodc/TreeLS')
#remotes::install_github('tiagodc/TreeLS')
#install.packages("TreeLS", dep = T)

packs_to_install <- c("lidR", "RANN", "conicfit", "alphashape3d", "alphahull",
                      "plotrix", "dbscan", "doParallel", "data.table",
                      "ks", "Metrics", "VoxR", "RCSF", "Morpho",
                      "ADPclust", "densityClust", "mgcv", "spatstat",
                      "flexclust", "matrixStats", "Distance", "lmfor",
                      "geosphere", "recexcavAAR", "DescTools", "foreach",
                      "benchmarkme", "stringi", "raster", "plyr", "png", "jpeg")



#install.packages("doMC", repos="http://R-Forge.R-project.org")
#devtools::install_github("omegahat/RDCOMClient")

install.packages(packs_to_install, dep = T)
# copy edci, Distance and mrds manually to R

# library("doParallel")
# library("data.table")
# library("ADPclust")
# library("densityClust")
# library("dae")
# library("plyr")
# library("doParallel")
# library("data.table")
# library("ADPclust")
# library("densityClust")
# library("dae")
# library("plyr")
# library("spatstat")
# library("alphahull")
# library("RANN")
# library("flexclust")
# library("sp")
# library("matrixStats")
# library("Distance")
# library("Rdistance")
# library("lmfor")
# library("rgl")
# library("conicfit")
# library("MASS")
# library("alphahull")
# library("igraph")
# library("geosphere")
# library("pracma")
# library("DescTools")
# library("mgcv")
# library("recexcavAAR")
# library("raster")
# library("rlas")
# library("lidR")
# library("TreeLS")
# library("dbscan")
#
