
.onLoad <- function(libname, pkgname) {
  # Access your package namespace
  ns <- asNamespace(pkgname)
  
  # Remove the 'knn' function if it exists
  if (exists("knn", envir = ns, inherits = FALSE)) {
    rm(list = "knn", envir = ns)
  }
}

library("plyr")
library("dplyr")
library("lidR")
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



.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    sprintf("%s version %s loaded.",
            pkgname,
            utils::packageVersion(pkgname))
  )
}

#' v.env stores the path _total_ground_veg
#' @export
v.env <- new.env(parent = emptyenv())
v.env$groundPath <- "_total_ground_veg/"
v.env$rndV <- "1989e10e48ae4bbc7a6ed615367367768ea2bd574929a0ea1fa118c031973596"
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

#' Separate function to generate the addressing string to store the files with all the settings in the name
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


#' Creates values with thousandmark (are prettier to look at)
#' @param val Numeric value that is extended with thousand mark 15200000 -> 15.200.000
#' @export
thMk <- function(val) {
  val2 <- format(val, big.mark = ".", decimal.mark = ",", scientific = FALSE)
  return(val2)
}




#' Convert abbreviation or numeric code to tree species
#' @param number abbreviation for tree species with two letters, be for beech, sp for spruce (English to German) 
#' @param modeNumber set to TRUE when the input is a number and not an abbreviation 
#' @export
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

#' Convert abbreviation or numeric code to tree species color
#' @param number abbreviation for tree species with two letters, be for beech, sp for spruce (English to German), or numeric
#' @param modeNumber set to TRUE when the input is a number and not an abbreviation 
#' @export
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
#' @param nr_plots_parallel how many plots should be calculated simultaneously
#' @param nr_cores_per_plot how many cores should each plot use
#' @param maximum_cores quick setting to process 20 plots at the same time with 10 cores per plot
#' @param detectTrees runs functions extractVegetation and clustSplit
#' @param segmentTrees runs functions crownFeel and computeTreeParams
#' @param crownParameters also performs stem analysis and extracts crown basal height, projection area and hull volume 
#' @param createAppFiles runs function createAppFiles with new background images
#' @param circleRadius draw additional circle in the background files
#' @param setting quick preset for 0=ground, 1=trees, 2=segment, a=app, 9=all
#' #' @export
processPlotsParallel <- function (inputFiles, fileFinders = "", 
                                  dirPath = paste0(getwd(), "/"),
                                  maximum_cores = 0,
                                  nr_plots_parallel = 20, 
                                  nr_cores_per_plot = 1,
                                  trafoFiles = "", 
                                  groundModels = TRUE, 
                                  detectTrees = TRUE, 
                                  segmentTrees = FALSE, 
                                  crownParameters = TRUE,
                                  createAppFiles = TRUE,
                                  setting = "",
                                  
                                  clip.trajectory.distance = 0, 
                                  clip.radius = 0, 
                                  
                                  pixelUnit_cm = 2, 
                                  circleRadius = 0, 
                                  drawTraj = F, 
                                  createBGR_pic = T, 
                                  
                                  voxelSize = 3, 
                                  limitShare = 0.003, 
                                  limitStems = 50,
                                  zScale = 3, 
                                  
                                  tileClipping = 2 # 2x2 for plots, 3x3 or 4x4 for large 1 ha scans
                                  
                                  )
{
  
  if(setting != "") cat("\n\nWARNING - SETTING NOT IMPLEMENTED YET, HAS NO EFFECT!!\n\n")
  library(digest)
  library(getPass)
  # error checking for number of core settings
  {
    nr_plots_parallel <- round(nr_plots_parallel)
    nr_cores_per_plot <- round(nr_cores_per_plot)
    if(nr_plots_parallel < 0){
      nr_plots_parallel <- 20
    }
    if(nr_cores_per_plot < 0){
      nr_cores_per_plot <- 1 # nested parallelization haut noch nicht so ganz hin...
    }
    if(maximum_cores == 1){
      nr_plots_parallel <- min(length(inputFiles), 200)
      nr_cores_per_plot <- 1
      # if(detectTrees == FALSE && segmentTrees == FALSE){
      #   nr_plots_parallel <- 200
      #   nr_cores_per_plot <- 1
      # }
    }
    cores_used <- nr_plots_parallel * nr_cores_per_plot
    if(detectTrees == FALSE && segmentTrees == FALSE){
      nr_plots_parallel <- cores_used
      nr_cores_per_plot <- 1
    }
    
    # PASSWORD PROTECTION FOR ALL CORES
    if(cores_used > 20){
      useAllCores <- F
      password <- getPass::getPass("Please confirm to access all WAFO cores: ", 
                                   forcemask = T)
      
      if (digest(password, algo = "sha256") == v.env$rndV) {
        cat("Access granted!\n")
        useAllCores <- T
      } else {
        password <- getPass::getPass("Wrong password, do you want to downgrade to 20 cores? (type yes)")
        if(is.null(password)){
          stop("User terminated script.")
        }
        if(is.element(password, c("yess", "yes", "y", "ye", "j", "ja", "z", "zes", "ze"))){
          useAllCores <- F
        } else if(digest(password, algo = "sha256") == v.env$rndV){
          cat("Access granted!\n")
          useAllCores <- T
        } else {
          stop("User terminated script.")
        }
      }
      if(!useAllCores){
        nr_plots_parallel <- 20
        nr_cores_per_plot <- 1
      }
    }
  }
  
 
  
  
  
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
    cat("#    ", nr_plots_parallel, " PLOTS PARALLEL\n", sep = "")
    cat("#    ", nr_cores_per_plot, " CORES PER PLOT\n", sep = "")
    cat("#  TODAY IS", paste(Sys.time()),"\n")
    cat("#  \n")
    
    
    if(groundModels){
      if(clip.trajectory.distance > 0){
        cat("#  trajectory clip at", clip.trajectory.distance, "m - ")
      }
      if(clip.radius > 0){
        cat("#  circle clip with radius of", clip.radius, "m - ")
      }
    }
    if(segmentTrees){
      cat("using", tileClipping*tileClipping, "tiles in", tileClipping, "x", tileClipping, "\n")
      cat("#  limitShare = ", limitShare, "  limitStems = ", limitStems, "%\n")
      cat("#      zScale = ", zScale, "x     voxelSize = ", voxelSize, " cm\n")
    }
    cat("#\n")
    
    trafoNums <- sum(trafoFiles == "")
    trafoNums <- length(trafoFiles) - trafoNums
    if(trafoNums == 0){
      trafoFiles <- rep("", length(inputFiles))
    } else {
      if(length(trafoFiles)!=length(inputFiles)){
        stop("Transformation files not provided for all input files.\nPlease double check trafo files and provide one for each plot!")
      }
      cat("# ", trafoNums, "GLOBAL TRANSFORMATION MATRICES PROVIDED\n")
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
  
  
  cat("Going parallel on", nr_plots_parallel, "x", nr_cores_per_plot, "cores for", length(fileFinders),"plots.\n")
  
  useProgressBar <- FALSE # works only for windows so far
  if (.Platform$OS.type == "windows") {
    library(doParallel)
    cl <- makeCluster(nr_plots_parallel)
    registerDoParallel(cl)
    cat("   (on windows system using doParallel)\n")
  } else {
    library(mcprogress)
    library(doMC)
    
    useProgressBar <- TRUE
    registerDoMC(cores = nr_plots_parallel)
    cat("   (on Unix-like system using doMC)\n")
    
  }
  cat("\n")
  
  
  library(foreach)
  start <- Sys.time()
  
  weightDetectTrees <- 30 # how much longer than compared to ground models
  weightSegmentTrees <- 8
  weightComputeTrees <- 13
  progressBarSteps <- sum(groundModels*1, detectTrees*weightDetectTrees, 
                          segmentTrees*weightSegmentTrees, 
                          crownParameters*weightComputeTrees, 
                          createAppFiles*1)
  consoleStartPath <- paste0(dirPath, "/parallel_console/", "par_", format(Sys.time(), "%y%m%d_%H%M"), "/")
  
  {
  foreach(i=1:length(fileFinders),  .errorhandling = 'remove', 
          .packages = c("treeX"), .combine = cbind) %dopar% {
            
            
            t1 <- Sys.time()
            nowStep <- 0
            if(useProgressBar){
              mcProgressBar(i, length(fileFinders), 
                            cores = nr_plots_parallel, 
                            subval = nowStep/progressBarSteps,
                            start = start)
            }
            nowLAZ <- inputFiles[i]
            fileFinder <- fileFinders[i]
            nowTrafo <- trafoFiles[i]
            
            #consolePath <- paste0(dirPath, "/parallel_console/", "par_", format(Sys.time(), "%y%m%d_%H%M"), "/")
            consolePath <- consoleStartPath
            if(!dir.exists(consolePath)) dir.create(consolePath, recursive = T)
            file_parallelProtocol <- paste0(consolePath, "/", fileFinder, "_par_Rcons.txt")
            file.create(file_parallelProtocol)
            sink(file_parallelProtocol, append = T)
            
            {
              #LASfile <- NA
              cat("\n#########################################\n",
                  "#\n",
                  "#  STARTING PARALLEL PROCESSING \n",
                  "#           FOR PLOT ", fileFinder, "\n",
                  "#           IN \"", getwd(), "\"\n", sep = "")
              cat("#  \n")
              cat("#    ", nr_plots_parallel, " PLOTS PARALLEL\n", sep = "")
              cat("#    ", nr_cores_per_plot, " CORES PER PLOT\n", sep = "")
              cat("#  TODAY IS", paste(Sys.time()),"\n")
              cat("#  \n")
              
              tasksScheduled <- paste0("#  ", 
                                       paste0(
                                         c(ifelse(groundModels, "GRD", "   "),
                                         ifelse(detectTrees, "DTC", "   "),
                                         ifelse(segmentTrees, "SEG", "   "),
                                         ifelse(crownParameters, "CRO", "   "),
                                         ifelse(createAppFiles, "APP", "   ")), 
                                         collapse = " + "
                                         ), 
                                       "\n")
              cat(tasksScheduled)
              cat("#  \n")
              
              
              if(groundModels){
                if(clip.trajectory.distance > 0){
                  cat("#  trajectory clip at", clip.trajectory.distance, "m - ")
                }
                if(clip.radius > 0){
                  cat("#  circle clip with radius of", clip.radius, "m - ")
                }
              } else {
                cat("#  \n")
              }
              if(segmentTrees){
                cat("using", tileClipping*tileClipping, "tiles in", tileClipping, "x", tileClipping, "\n")
                cat("#  limitShare = ", limitShare, "  limitStems = ", limitStems, "%\n")
                cat("#      zScale = ", zScale, "x     voxelSize = ", voxelSize, " cm\n")
              }
              cat("#\n")
              
              if(nowTrafo != ""){
                cat("#  GLOBAL TRANSFORMATION MATRICE PROVIDED\n")
                cat("#\n")
              }
              cat("#########################################\n")
              cat("\n\n")
            }
            
            
            if(groundModels){
              try(extractVegetation(nowLAZ, fileFinder, 
                                    trafoMatrix.path = nowTrafo, 
                                    clip.trajectory.distance = clip.trajectory.distance, 
                                    clip.radius = clip.radius))
              if(useProgressBar){
                nowStep <- nowStep+1
                mcProgressBar(i, length(fileFinders), 
                              cores = nr_plots_parallel, 
                              subval = nowStep/progressBarSteps,
                              start = start)
              }
            }
            
            if(detectTrees){
              try(clustSplit(fileFinder = fileFinder, filterINT = 97, 
                             nr_cores = nr_cores_per_plot, 
                             retainPointClouds = T))
              if(useProgressBar){
                nowStep <- nowStep+weightDetectTrees
                mcProgressBar(i, length(fileFinders), 
                              cores = nr_plots_parallel, 
                              subval = nowStep/progressBarSteps,
                              start = start)
              }
            }
            
            if(segmentTrees){
              try(crownFeel(fileFinder, 
                            voxelSize = voxelSize, 
                            limitShare = limitShare, 
                            limitStems = limitStems, 
                            zScale = zScale,
                            tileClipping = tileClipping, frame.rad = 1.1, 
                            retainPointClouds = T))
              if(useProgressBar){
                nowStep <- nowStep+weightSegmentTrees
                mcProgressBar(i, length(fileFinders), 
                              cores = nr_plots_parallel, 
                              subval = nowStep/progressBarSteps,
                              start = start)
              }
              try(computeTreeParams(fileFinder, getRAM = F, 
                                    voxelSize = voxelSize,
                                    limitShare = limitShare, 
                                    zScale = zScale, nr_cores = 1, nr_cores_params = nr_cores_per_plot,
                                    writeLAZ = T, 
                                    retainPointClouds = F, 
                                    crownParameters = crownParameters))
              if(useProgressBar){
                nowStep <- nowStep+weightComputeTrees
                mcProgressBar(i, length(fileFinders), 
                              cores = nr_plots_parallel, 
                              subval = nowStep/progressBarSteps,
                              start = start)
              }
            } else if(crownParameters){
              try(computeCrownParams(fileFinder, 
                                    voxelSize = voxelSize,
                                    limitShare = limitShare, 
                                    zScale = zScale, nr_cores = nr_cores_per_plot, 
                                    retainPointClouds = F))
            }
            
            if(createAppFiles){
              try(createAppFiles(fileFinder,
                                 pixelUnit_cm = pixelUnit_cm, eraseSpecies = T, drawTraj = drawTraj, 
                                 #drawLines = drawLines, slices = nowSlices,
                                 #isoLines = 1,
                                 createBGR_pic = createBGR_pic, 
                                 circleRadius = circleRadius, 
                                 drawGround = T, 
                                 drawRedSlice = T))
              if(useProgressBar){
                nowStep <- nowStep+1
                mcProgressBar(i, length(fileFinders), 
                              cores = nr_plots_parallel, 
                              subval = nowStep/progressBarSteps,
                              start = start)
              }
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
  if(useProgressBar) closeProgress(start)
  }
}





#' 
#' Copies (or moves) all laz and trajectory files
#' from the FARO connect directory to a specified path
#' 
#' 
#' @param inConnectPath path to the input folders of all scans
#' @param outNewFolder path where all laz and trajs will be copied (or moved)
#' @param move set TRUE to move instead of copy files
#' @param handleCleanLAZ what to do with files ending in _clean.laz: -1 (default) ignore them, 0 move them as well, 1 only move _clean.laz files
#' @param cleanConnectPath set TRUE to delete all other files in inConnectPath
#' @export
pullFromFaroConnect <- function(inConnectPath, outNewFolder, 
                                move = F, 
                                handleCleanLAZ = -1, 
                                cleanConnectPath = F){
  
  
  # Hier Ordner mit FARO Connect Dateien
  inConnectPath <- "D:/Zoran/_in_connect/"
  
  # Wohin verschieben?
  outNewFolder <- "D:/Zoran/_in_laz/"
  
  if(file.exists(paste0(inConnectPath, "g/sheet/tsc/"))){
    inConnectPath <- paste0(inConnectPath, "g/sheet/tsc/")
  }
  cFiles <- list.files(inConnectPath, recursive = T)
  cFiles <- cFiles[endsWith(cFiles, "gs-traj") | endsWith(cFiles, "laz")]
  if(handleCleanLAZ <= -1){
    cFiles <- cFiles[!endsWith(cFiles, "_clean.laz")]
  } else if(handleCleanLAZ >= 1){
    cFiles <- cFiles[!duplicated(gsub("_clean.laz", ".laz", cFiles), fromLast = T)]
  }
  {
    if(move){
      cat("Moving ")
    } else {
      cat("Copying ")
    }
    cat(length(cFiles), "files\n from: ", inConnectPath, "\n   to: ", outNewFolder)
  }
  rento <- gsub("ts/", "", basename(cFiles))
  rento <- gsub("ts../", "", rento)
  rento <- gsub("ts.../", "", rento)
  rento
  
  # Loop through each subdirectory
  for (i in 1:length(cFiles)) {
    old_dir <- paste0(inConnectPath, cFiles[i])
    
    # Define the new directory path in the target directory
    new_dir <- paste0(outNewFolder, "", rento[i])
    
    cat(i, "- file", old_dir, "\n")
    #if(endsWith(new_dir, "geoslam")) next()
    
    #if(!endsWith(new_dir, "gs-traj")) next()
    
    if(!file.exists(old_dir)) next()
    
    if(!dir.exists(dirname(new_dir))) dir.create(dirname(new_dir), recursive = T)
    
    # Get the original modification time of the subdirectory
    orig_time <- file.info(old_dir)$mtime
    
    # Copy or Move the subdirectory to the target directory
    file.copy(from = old_dir, to = new_dir, overwrite = T)
    # MOVE - careful, old file is lost
    # file.renamefrom = old_dir, to = new_dir, overwrite = T)
    
    #Sys.setFileTime(dirname(old_dir), orig_time)
    # Set the modification time for the new directory
    if (length(orig_time) > 0) {
      Sys.setFileTime(new_dir, orig_time)
      Sys.setFileTime("D:/Laser_Tirol_Noflatscher/_in_laz/", orig_time)
    }
    
    # Deleting processed .geoslam file and rest of folder
    if(cleanConnectPath) unlink(paste0(old_dir, "/"), recursive = T)
    
  }
}
