library("lidR")
library("TreeLS")
library("RANN")
library("conicfit")
library("alphashape3d")
library("alphahull")
library("plotrix")
library("dbscan")

if (.Platform$OS.type == "windows") {
  library(doParallel)
} else {
  library(doMC)
}

library("plyr")
library("sf")
library("sp")
library("spatstat")
library("Morpho") # for transformation of point clouds




v.env <- new.env(parent = emptyenv())
v.env$groundPath <- "_total_ground_veg/"
# groundPath <- v.env$groundPath


#groundPath <<- "_total_ground_veg/"


# if(!exists("dirPath")){
#   dirPath <- paste0(getwd(),"/")
#   # sets the global path for all output files (vegetation, trees, lists)
# }
#dirPath lost its meaning, as I can't access the file within the package,
#need to accomodate with the thought of using setwd() and getwd()


# dirPath <- paste0(getwd(),"/")
# cat("Working directory (dirPath) for all files is:",dirPath,"\n")
# #cat("Change it in the variable dirPath.\n")
# if(!dir.exists(paste0(dirPath,groundPath))) dir.create(paste0(dirPath,groundPath))
# memory.size(TRUE) #



# GLOBAL SETTINGS (can be modified to fit best to every set)
if(1==2){
  # only for Zenodo files at D: hard drive
  chooseFile <- function(number){
    if(number > 0 && number < 20){
      readInput <<- paste0("D:/zenodo/PLS_clouds/",number,".laz")     # full path and name of the laser point cloud in dirPath
      fileFinder <<- paste0("p",number)
      try(ref.plot_id <<- strtoi(substr(fileFinder,2,999)))

      cat("Choosing active point cloud",fileFinder,"\n")
      if(!file.exists(readInput)){
        cat("ERROR - file",readInput,"not found!\n")
        readInput <<- NULL
        fileFinder <<- NULL
      }
    }
  }

  #bushPreparation <- FALSE
  #filterSOR <- FALSE
  if(!exists("groundCutHeight")) groundCutHeight <<- 0.50 #50 cm of the stems are lost due to ground removal (DTM model)

  if(!exists("ref")) ref <<- "D:/zenodo/meta/reference_data.csv"

  if(!exists("readInput")) readInput <<- paste0("D:/zenodo/PLS_clouds/4.laz")     # full path and name of the laser point cloud in dirPath
  if(!exists("fileFinder")) fileFinder <<- "p4"
  if(!exists("ref.plot_id")) ref.plot_id <<- 4

  if(!exists("mode")) mode <- ""
  if(!exists("clipHeight")) clipHeight <<- 3
  if(!exists("bottomCut")) bottomCut <<- 1
  if(!exists("numberOfPoints")) numberOfPoints <<- 300
  if(!exists("heightExtent")) heightExtent <<- 1.3

  if(!exists("cutWindow")) cutWindow <<- c(-1000,-1000,2000)
  # cutWindow <- c(-2,-4,5)
  # cutWindow <- c(2,-15,5)

  if(!exists("filterSOR")) filterSOR <<- FALSE

  if(!exists("zScale")) zScale <<- 2 # double the height progress
  if(!exists("limitShare")) limitShare <<- 0.004
  if(!exists("totalRuns")) totalRuns <<- 1000

  if(!exists("silent")) silent <<- TRUE
  if(!exists("fast")) fast <<- FALSE
  if(!exists("retainPointClouds")) retainPointClouds <<- TRUE
}



# FUNCTIONS GLOBALLY

# separate function to generate the adressing string to store the files with all the settings in the name
#' @export
generateSetString <- function(fileFinder = NA, mode = NA, threshold = NA, threshold_percent = NA,
                              bottomCut = NA, clipHeight = 0, bushPreparation = FALSE,
                              filterSOR = FALSE, level = NA, numberOfPoints = NA,
                              cutWindow = c(-1000,-1000,2000),
                              limitShare = NA, zScale = NA, voxelSize = 0,
                              slices = NA, bottom = NA, step = NA, overlap = NA,
                              referenced = FALSE,
                              silent = FALSE){
  # MODES: way of finding the initial tree centers
  #        ___COMP (old approach Cloud Compare Connected Components)
  #        ___SLIS (slice sampling BSC thesis)
  setString <- ""
  if(!is.na(fileFinder)){
    setString <- paste0(fileFinder)
    if(!silent){
      cat("These are your settings:\n", sep = "")
      cat("     We are located on ",fileFinder,"\n", sep = "")
    }
  }
  if(!is.na(mode)){
    if(!silent)cat("     the mode will be",mode,"\n")
    setString <- paste0(setString,"_",mode)
  }
  if(!is.na(threshold_percent)){
    if(!silent)cat("     using only the",threshold_percent,"% of maximum intensity values\n")
    setString <- paste0(setString,"_",threshold_percent)
  } else if(!is.na(threshold)){
    if(!silent)cat("     using only the",thMk(threshold),"maximum intensity values\n")
    setString <- paste0(setString,"_",threshold,"pts")
  }

  if(!is.na(slices)){
    setString <- paste0(setString,"_x",slices)
  }
  if(!is.na(bottom)){
    setString <- paste0(setString,"_",bottom,"cm")
  }
  if(!is.na(step)){
    setString <- paste0(setString,"_inc",step)
  }
  if(!is.na(overlap)){
    setString <- paste0(setString,"+",overlap)
  }

  if(clipHeight > 0){
    setString <- paste0(setString, "_",round(bottomCut*100),"to",round(clipHeight*100))
  }
  if(bushPreparation){
    setString <- paste0(setString, "b")
  }
  if(filterSOR){
    #SOR NOT IMPLEMENTED YET - FIND AN R SOLUTION FOR IT IF NECCESSARY
    #setString <- paste0(setString, "SOR")
  }
  if(!is.na(level)) setString <- paste0(setString,"_l",level)
  if(!is.na(numberOfPoints)) setString <- paste0(setString,"a",numberOfPoints)


  if(!is.na(limitShare)) setString <- paste0(setString,"_lim",limitShare)
  if(!is.na(zScale)) setString <- paste0(setString, "_Z",zScale)
  if(voxelSize != 0) setString <- paste0(setString, "_vox",voxelSize)

  if(sum(cutWindow == c(-1000,-1000,2000))!=3){
    try({
      XL <- cutWindow[1]
      YL <- cutWindow[2]
      width <- cutWindow[3]
      locationStr <- paste0("x",XL,"y",YL,"dim",width)

      if(setString != ""){
        setString <- paste0(setString, "_",locationStr)
      } else {
        setString <- locationStr
      }
    })
  }
  if(referenced) setString <- paste0(setString, "_ref")


  return(setString)
}


#Creates values with thousandmark (are prettier to look at)
#' @export
thMk <- function(val) {
  val2 <- format(val, big.mark = ".", decimal.mark = ",", scientific = FALSE)
  return(val2)
}



treeSpecies <- function(number, modeNumber = FALSE) {
  species <- "unklar"
  if(modeNumber){
    species <- switch (paste0(number),
                       "10" = "Buche",
                       "1" = "Fichte",
                       "2" = "Tanne",
                       "3" = "Laerche",
                       "11" = "Eiche",
                       "4" = "Kiefer",
                       "25" = "Pappel",
                       "27" = "Weide",
                       "19" = "Birke",
                       "12" = "Hainbuche",
                       "18" = "Kirsche",
                       "14" = "Ahorn",
                       "15" = "Ulme",
                       "13" = "Esche",
                       "unbekannt")
  } else {
    species <- switch (paste0(number),
                       "be" = "Buche",
                       "sp" = "Fichte",
                       "fir" = "Tanne",
                       "la" = "Laerche",
                       "oak" = "Eiche",
                       "pin" = "Kiefer",
                       "pib" = "Schwarzkiefer",
                       "con" = "sonst. NH",
                       "pop" = "Pappel",
                       "sall" = "Weide",
                       "bir" = "Birke",
                       "hb" = "Hainbuche",
                       "cher" = "Kirsche",
                       "map" = "Ahorn",
                       "sor" = "Sorbus",
                       "elm" = "Ulme",
                       "ash" = "Esche",
                       "unbekannt")
  }
  return(species)
}

treeCol <- function(number, modeNumber = FALSE) {
  color <- "pink"
  if(modeNumber){
    color <- switch (paste0(number),
                     "1" = "chartreuse",
                     "2" = "darkgreen",
                     "3" = "darkgoldenrod1",
                     "4" = "cyan",
                     "10" = "blue",
                     "11" = "darkred",
                     "25" = "brown",
                     "27" = "cadetblue1",
                     "19" = "aquamarine",
                     "12" = "cornsilk3",
                     "18" = "deeppink4",
                     "14" = "darkorange1",
                     "15" = "yellow",
                     "13" = "burlywood4",
                     "pink")
  } else {
    color <- switch (paste0(number),
                     "sp" = "chartreuse",
                     "fir" = "darkgreen",
                     "la" = "darkgoldenrod1",
                     "pin" = "cyan",
                     "pib" = "darkolivegreen",
                     "con" = "darkseagreen",
                     "be" = "blue",
                     "oak" = "darkred",
                     "pop" = "brown",
                     "sall" = "cadetblue1",
                     "bir" = "aquamarine",
                     "hb" = "cornsilk3",
                     "cher" = "deeppink4",
                     "map" = "darkorange1",
                     "sor" = "brown1",
                     "elm" = "yellow",
                     "ash" = "burlywood4",
                     "pink")
  }
  return(color)
}






#' @export
centrit <- function(pts){
  centre = st_coordinates(st_centroid(st_union(pts)))
  pts = st_coordinates(pts)

  theta = atan2(pts[,2]-centre[,2], pts[,1]-centre[,1])
  poly = cbind(
    pts[order(theta),1],
    pts[order(theta),2])
  poly = st_polygon(list(rbind(poly, poly[1,])))
  poly

}

#' @export
bestMclust <- function (clust, nc = 1, crit = "value")
{
  c = class(clust)
  if (!is.matrix(clust)) {
    result = clust
  }
  else {
    if (!(crit %in% colnames(clust))) {
      #cat("No row \"", crit, "\" found!\n")
      result = clust
    }
    else {
      result = rbind(clust[order(clust[, crit], decreasing = TRUE)[1:min(nc,
                                                                         nrow(clust))], ])
      class(result) = c
    }
  }
  result
}

# clust = cbind(cx, cy, r, value)[1:count, ]
# prec = prec
# ncol = 3
# dz = TRUE
#' @export
deldupMclust <- function (clust, prec = NULL, ncol = NULL, dz = TRUE)
{
  c = class(clust)
  if(length(c) > 1){
    c <- c[1]
  }
  if (is.matrix(clust) && dz && any(colnames(clust) == "value")) {
    if (sum(clust[1:nrow(clust), "value"] != 0) == 0)
      clust = clust[1, ]
    else clust = clust[1:nrow(clust) * (clust[1:nrow(clust),
                                              "value"] != 0), ]
  }
  if (is.matrix(clust)) {
    if (all(is.null(dimnames(clust))) | all(dimnames(clust)[[2]] !=
                                            "count")) {
      count = rep(1, nrow(clust))
      clust = cbind(clust, count)
    }
    if (is.null(ncol)) {
      if (!is.null(c) && c == "oregMclust")
        ncol = 2
      else if (!is.null(c) && c == "circMclust")
        ncol = 3
      else ncol = ncol(clust) - sum(colnames(clust) %in%
                                      c("value", "count", "proj", "env"))
    }
    if (ncol > ncol(clust))
      ncol = ncol(clust)
    if (!is.null(prec)) {
      clust[, 1:ncol] = round(clust[, 1:ncol], prec)
      if (!is.null(c) && c == "oregMclust")
        clust[, 1] = clust[, 1]%%round(2 * pi, prec)
    }
    result = clust
    count = 1
    for (i in 2:nrow(clust)) {
      j = 1
      ins = TRUE
      while (j <= count && ins) {
        if (is.na(all(clust[i, 1:ncol] == result[j,
                                                 1:ncol])) | all(clust[i, 1:ncol] == result[j,
                                                                                            1:ncol])) {
          ins = FALSE
          if (!is.na(all(clust[i, 1:ncol] == result[j,
                                                    1:ncol])))
            result[j, "count"] = result[j, "count"] +
              clust[i, "count"]
        }
        j = j + 1
      }
      if (ins) {
        count = count + 1
        result[count, ] = clust[i, ]
      }
    }
    result = result[1:count, ]
    if (count > 1)
      class(result) = c
  }  else result = clust
  result
}






# prec = 4
# minsx = min(datax)
# maxsx = max(datax)
# minsy = min(datay)
# maxsy = max(datay)
# minsr = 0.01 * max(datax, datay)
# maxsr = (max(datax, datay) - min(datax, datay))
# nsc = 5
# nc = NULL
#
# minsd = NULL
# maxsd = NULL
# brminx = minsx
# brmaxx = maxsx
#
# brminy = minsy
# brmaxy = maxsy
# brminr = minsr
# brmaxr = maxsr
#
# brmaxit = 1000
#
# datax=coords[,1]
# datay=coords[,2]
# method="const"
# # nx=1, ny=1, nr=1, # nx=10, ny=10, nr=25,
# nx=10
# ny=10
# nr=5 #nx=25, ny=25, nr=5,
# minsr=0.01
# maxsr=0.5
# nc=1 # nc=1,
# # minsd=0.1, maxsd=0.5,
# bw=0.01 # bw=0.05
#


#' @export
circMclust <- function (datax, datay, bw, method = "const", prec = 4, minsx = min(datax),
                        maxsx = max(datax), nx = 10, minsy = min(datay), maxsy = max(datay),
                        ny = 10, minsr = 0.01 * max(datax, datay), maxsr = (max(datax,
                                                                                datay) - min(datax, datay)), nr = 10, nsc = 5, nc = NULL,
                        minsd = NULL, maxsd = NULL, brminx = minsx, brmaxx = maxsx,
                        brminy = minsy, brmaxy = maxsy, brminr = minsr, brmaxr = maxsr,
                        brmaxit = 1000)
{
  #cat("Break with <CTRL>-C (linux) or <ESC> (windows)\n")
  n = min(length(datax), length(datay))
  count = integer(1)
  if (method == "all") {
    m = 0
    nmax = n * (n - 1)/2
    nmax = nmax * (nmax - 1)/2
    cx = double(nmax)
    cy = double(nmax)
    r = double(nmax)
    value = double(nmax)
    result = .C("c_oregMcirc", as.double(datax), as.double(datay),
                n, as.double(bw), as.integer(m), 0L, 0L, 0L, as.double(maxsr),
                0, 0, 0, 0L, as.double(brminx), as.double(brmaxx),
                as.double(brminy), as.double(brmaxy), as.double(brminr),
                as.double(brmaxr), as.integer(brmaxit), cx = cx,
                cy = cy, r = r, value = value, count = count)
    cx = result$cx
    cy = result$cy
    r = result$r
    value = result$value
    count = result$count
  } else if (method == "prob") {
    m = 1
    nmax = 100
    cx = double(nmax)
    cy = double(nmax)
    r = double(nmax)
    value = double(nmax)
    nocirc = -1
    count = 0
    cxnew = double(1)
    cynew = double(1)
    rnew = double(1)
    vnew = double(1)
    c = integer(1)
    if (!is.null(minsd))
      minsdq = minsd^2
    if (!is.null(maxsd))
      maxsdq = maxsd^2
    repeat {
      i = 1
      while (i <= nsc) {
        repeat {
          i1 = round(runif(1, min = 1, max = n))
          i2 = round(runif(1, min = 1, max = n))
          i3 = round(runif(1, min = 1, max = n))
          if (i1 != i2 && i1 != i3 && i2 != i3)
            break
        }
        if ((is.null(minsd) || (((datax[i1] - datax[i2])^2 +
                                 (datay[i1] - datay[i2])^2) > minsdq && ((datax[i1] -
                                                                          datax[i3])^2 + (datay[i1] - datay[i3])^2) >
                                minsdq && ((datax[i3] - datax[i2])^2 + (datay[i3] -
                                                                        datay[i2])^2) > minsdq)) && (is.null(maxsd) ||
                                                                                                     (((datax[i1] - datax[i2])^2 + (datay[i1] -
                                                                                                                                    datay[i2])^2) < maxsdq && ((datax[i1] -
                                                                                                                                                                datax[i3])^2 + (datay[i1] - datay[i3])^2) <
                                                                                                      maxsdq && ((datax[i3] - datax[i2])^2 + (datay[i3] -
                                                                                                                                              datay[i2])^2) < maxsdq))) {
          result = .C("c_oregMcirc", as.double(datax),
                      as.double(datay), n, as.double(bw), as.integer(m),
                      as.integer(i1 - 1), as.integer(i2 - 1),
                      as.integer(i3 - 1), as.double(maxsr), 0,
                      0, 0, 0L, as.double(brminx), as.double(brmaxx),
                      as.double(brminy), as.double(brmaxy), as.double(brminr),
                      as.double(brmaxr), as.integer(brmaxit),
                      cxnew = cxnew, cynew = cynew, rnew = rnew,
                      vnew = vnew, c = c)
          cxnew = result$cxnew
          cynew = result$cynew
          rnew = result$rnew
          vnew = result$vnew
          c = result$c
          if (c > 0) {
            count = count + 1
            if (length(cx) < count) {
              cx = c(cx, double(nmax))
              cy = c(cy, double(nmax))
              r = c(r, double(nmax))
              value = c(value, double(nmax))
            }
            cx[count] = cxnew
            cy[count] = cynew
            r[count] = rnew
            value[count] = vnew
          }
          i = i + 1
        }
      }
      nocircnew = nrow(deldupMclust(cbind(cx, cy, r)[1:count,
      ], prec = prec))
      if (is.null(nocircnew))
        nocircnew = 1
      #cat("Found clusters: ", nocircnew, "\n")
      if (is.null(nc) && nocirc == nocircnew)
        break
      if (!is.null(nc))
        if (nocircnew >= nc)
          break
      nocirc = nocircnew
    }
  } else if (method == "const") {
    m = 2
    nmax = 100
    cx = double(nmax)
    cy = double(nmax)
    r = double(nmax)
    value = double(nmax)
    count = integer(1)
    circ = NULL
    xseq = seq(minsx, maxsx, length = nx)
    yseq = seq(minsy, maxsy, length = ny)
    rseq = seq(minsr, maxsr, length = nr)
    startpoints = matrix(c(rep(xseq, each = ny * nr), rep(rep(yseq,
                                                              each = nr), nx), rep(rseq, nx * ny)), ncol = 3)
    for (i in 1:(ceiling(nrow(startpoints)/nmax))) {
      s = seq(((i - 1) * nmax + 1), (min(i * nmax, nrow(startpoints))))
      # cat("calculating clusters for startingvalues ",
      #     s[1], "-", s[length(s)], " of ", nrow(startpoints),
      #     "...\n")
      result = .C("c_oregMcirc", as.double(datax), as.double(datay),
                  n, as.double(bw), as.integer(m), 0L, 0L, 0L,
                  as.double(maxsr), as.double(startpoints[s, 1]),
                  as.double(startpoints[s, 2]), as.double(startpoints[s,
                                                                      3]), as.integer(length(s)), as.double(brminx),
                  as.double(brmaxx), as.double(brminy), as.double(brmaxy),
                  as.double(brminr), as.double(brmaxr), as.integer(brmaxit),
                  cx = cx, cy = cy, r = r, value = value, count = count)
      cx = result$cx
      cy = result$cy
      r = result$r
      value = result$value
      count = result$count
      if (count > 0) {
        if (count == 1)
          circnew = c(cx[1], cy[1], r[1], value[1],
                      count)
        else circnew = deldupMclust(cbind(cx, cy, r,
                                          value)[1:count, ], prec = prec, ncol = 3)
        if (i == 1)
          circ = circnew
        else circ = deldupMclust(rbind(circ, circnew),
                                 prec = prec, ncol = 3)
      }
    }
    if (is.null(circ)) {
      #cat("no cluster found!\n")
    }
    else {
      #cat("finished\n")
      if (is.matrix(circ))
        circ[, 4] = -circ[, 4]
      else circ[4] = -circ[4]
    }
  }   else {
    #cat("unknown method\n")
  }
  if (method != "const") {
    value = -value
    circ = cbind(cx, cy, r, value)[1:count, ]
    circ = deldupMclust(circ, prec = prec, ncol = 3)
  }
  if (!is.null(circ)) {
    class(circ) = "circMclust"
    rownames(circ) = NULL
  }
  circ
}



#' Full tree detection and instance segmentation workflow
#'
#' Runs the complete workflow for multiple plots parallel. 
#' Be careful, as this process can be computationally demanding
#' and the system might fail when limited in RAM. 
#' I recommend having at least 20 GB of RAM for each sample plot radius 20 m
#' and at least 60 GB of RAM for each scan up to 1 ha. 
#'
#' @param inputFiles vector of paths to .laz or .las files as input for processing (plots)
#' @param fileFinders vector of file names for renaming input files. Must be of same length as lazFiles
#' @param trafoFiles vector of path to transformation matrices for each plot
#' @param nr_cores_plots how many plots should be calculated simultaneously
#' @param detectTrees runs functions extractVegetation and clustSplit
#' @param segmentTrees runs functions crownFeel and computeTreeParams
#' @param crownParameters also performs stem analysis and extracts crown basal height, projection area and hull volume 
#' @param createAppFiles runs function createAppFiles with new background images
#' #' @export
processPlotsParallel <- function (inputFiles, fileFinders = "", 
                                  dirPath = paste0(getwd(), "/"),
                                  nr_cores_plots = 4, trafoFiles = "", 
                                  detectTrees = TRUE, 
                                  segmentTrees = FALSE, 
                                  crownParameters = TRUE,
                                  createAppFiles = TRUE,
                                  
                                  clip.trajectory.distance = 0, 
                                  clip.radius = 0, 
                                  
                                  
                                  voxelSize = 3, 
                                  limitShare = 0.003, 
                                  limitStems = 50,
                                  zScale = 3, 
                                  
                                  tileClipping = 2 # 2x2 for plots, 3x3 or 4x4 for large 1 ha scans
                                  
                                  )
{
  
  dirPath <- paste0(dirPath, "/")
  dirPath <- gsub("//", "/", dirPath)
  if(!dir.exists(dirPath)) dir.create(dirPath)
  setwd(dirPath) # all the files and folder structure will be created here
  
  
  
  {
    #LASfile <- NA
    cat("\n#########################################\n",
        "#\n",
        "#  STARTING PARALLEL PROCESSING \n",
        "#           FOR ", length(inputFiles), " PLOTS\n",
        "#           IN \"", getwd(), "\"\n", sep = "")
    cat("#  \n")
    cat("#  TODAY IS", paste(Sys.time()),"\n")
    cat("#  \n")
    
    
    if(detectTrees){
      if(clip.trajectory.distance > 0){
        cat("#  clip trajectory to ", clip.trajectory.distance, "m\n")
      }
      if(clip.radius > 0){
        cat("#  clip circle of", clip.radius, "m\n")
      }
      cat("#  \n")
    }
    if(segmentTrees){
      cat("#   voxelSize = ", voxelSize, " cm\n")
      cat("#  limitShare = ", limitShare, "\n")
      cat("#  limitStems = ", limitStems, "%\n")
      cat("#      zScale = ", zScale, "x\n") 
      cat("#    and", tileClipping*tileClipping, "tiles are done", tileClipping, "x", tileClipping, "\n")
    }
    cat("#\n")
    
    if(trafoFiles[1] == ""){
      trafoFiles <- rep("", length(inputFiles))
    } else {
      if(length(trafoFiles)!=length(inputFiles)){
        stop("Transformation files not provided for all input files.\nPlease double check trafo files and provide one for each plot!")
      }
      cat("#  GLOBAL TRANSFORMATION MATRICES PROVIDED\n")
      cat("#\n")
    }
    cat("#########################################\n")
    cat("\n\n")
  }
  
  
  if(fileFinders[1] == ""){
    cat("Extracting filenames from input files:\n")
    fileFinders <- gsub(".laz", "", basename(inputFiles))
    fileFinders <- gsub(".las", "", fileFinders)
    cat(paste0(sprintf(" %03d", c(1:length(fileFinders))), "  ", fileFinders, collapse = "\n"))
    cat("\n\n")
  }
  
  
  
  cat("Going parallel on", nr_cores_plots, "cores.\n")
  
  if (.Platform$OS.type == "windows") {
    library(doParallel)
    cl <- makeCluster(nr_cores_plots)
    registerDoParallel(cl)
    cat("   (on windows system using doParallel)\n")
  } else {
    library(doMC)
    registerDoMC(cores = nr_cores_plots)
    cat("   (on Unix-like system using doMC)\n")
  }
  cat("\n")
  
  library(foreach)
  
  foreach(i=1:length(fileFinders),  .errorhandling = 'remove', 
          .packages = c("treeX"))%dopar% {
            
            t1 <- Sys.time()
            
            nowLAZ <- inputFiles[i]
            fileFinder <- fileFinders[i]
            trafoFiles <- trafoFiles[i]
            
            consolePath <- paste0(dirPath, "/parallel_console/")
            if(!dir.exists(consolePath)) dir.create(consolePath)
            file_parallelProtocol <- paste0(consolePath, fileFinder, "_par.txt")
            file.create(file_parallelProtocol)
            sink(file_parallelProtocol, append = T)
            
            
            if(detectTrees){
              
              try(extractVegetation(nowLAZ, fileFinder, 
                                    trafoMatrix.path = trafoFiles, 
                                    clip.trajectory.distance = clip.trajectory.distance, 
                                    clip.radius = clip.radius))
              
              try(clustSplit(fileFinder = fileFinder, filterINT = 97, 
                             nr_cores = 1, 
                             retainPointClouds = T))
              
              
            }
            
            if(segmentTrees){
              
              try(crownFeel(fileFinder, 
                            voxelSize = voxelSize, 
                            limitShare = limitShare, 
                            limitStems = limitStems, 
                            zScale = zScale,
                            tileClipping = tileClipping, frame.rad = 1.1, 
                            retainPointClouds = T))
              
              try(computeTreeParams(fileFinder, getRAM = F, 
                                    voxelSize = voxelSize,
                                    limitShare = limitShare, 
                                    zScale = zScale, nr_cores = 1, nr_cores_params = 1,
                                    writeLAZ = T, 
                                    retainPointClouds = F, 
                                    crownParameters = crownParameters))
              
            }
            
            if(createAppFiles){
              try(createAppFiles(fileFinder,
                                 pixelUnit_cm = 2, eraseSpecies = T, drawTraj = F, 
                                 #drawLines = drawLines, slices = nowSlices,
                                 #isoLines = 1,
                                 writeColoredLAZ = F, createBGR_pic = T, createJPG = T,
                                 drawGround = T, 
                                 drawRedSlice = T))
            }
            
            {
              t2 <- Sys.time()
              dt <- difftime(t2,t1)
              dt <- paste(round(dt, 1), units(dt))
              
              cat("#################################\n")
              cat("#\n")
              cat("#  Complete Workflow done for", fileFinder,"\n")
              cat("#  Time required: ", dt,"\n")
              cat("#\n")
              cat("#################################\n")
              cat("\n\n")
            }
            sink()
            
          }
  
  
}
