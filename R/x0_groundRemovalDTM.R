
if(!exists("LAS_veg")){
  LAS_veg <<- NA
  LAS_ground <<- NA
  LAS_veg_name <<- "blank"
  }

# no longer supported in later R-versions
# memory.limit(size = 100000)




#' Ground classification of raw point cloud with transformation into old set
#' for more help see function extractVegetation, 
#' which is run twice (local & transformed) within this function
#'
#' Reads in a .las or .laz file from disk and extracts the vegetation cloud
#' using the cloth simulation filter (scf) via lasground from the lidR package
#' automatically computing transformation matrix 
#' to merge new set fileFinder into old set trafo.matchOldSet
#'
#' @param trafo.matchDirPath directory of the old fileFinder which is used for
#' @param trafo.matchOldSet name of set to try for matching new set into
#' @param trafo.matchSlice can be set to cut a circle from 0 0 with the given radius of the big input file
#' @export
transformVegetation <- function(LASfile, fileFinder, 
                                match.voxel.new = 0, 
                                match.voxel.old = 0, 
                                match.voxel.both = 0,
                                runInitialExtractVegetation = T, 
                                groundMergeCut = 0, ipad = FALSE,
                                groundCutHeight = 1.0, steepSlope = TRUE, clothSize = 0.10,
                                doFilter = FALSE, filterRes = 0.05, filterN = 27, 
                                clip.radius = 0, clip.trajectory.distance = 0,
                                clip.x = 0, clip.y = 0,
                                draw.trajectory = TRUE, 
                                selector = "xyzcit0RGB", dtm.path = "",
                                preTrafoMatrix.path = "",
                                trafo.matchOldSet = "", 
                                trafo.matchDirPath = "",
                                trafo.matchSlice = "_clusterSlice_120to140.laz",
                                rainFilter = 0,
                                additionalSlices = TRUE,
                                exportSlice.upperLimit = 3, exportSlice.lowerLimit = 1, #m
                                exportClippedLAS = FALSE,
                                dirPath = paste0(getwd(), "/")){
  
  library(Morpho)
  
  if(match.voxel.both > 0){
    match.voxel.new <- match.voxel.both
    match.voxel.old <- match.voxel.both
  }
  if(trafo.matchOldSet != ""){
    if(trafo.matchDirPath == "") trafo.matchDirPath <- dirPath
    
    path_oldSlice <- paste0(trafo.matchDirPath, "_total_ground_veg/", 
                            trafo.matchOldSet, trafo.matchSlice)
    
    if(!file.exists(path_oldSlice)){
      stop("Can't transform without file ", path_oldSlice, "!")
    } 
    cat("TRANSFORMING SET", fileFinder, "\n  into set", trafo.matchOldSet, "\n  via file", path_oldSlice, "\n\n\n")
  }
  
  
  if(runInitialExtractVegetation){
    extractVegetation(LASfile = LASfile, fileFinder = fileFinder, 
                      groundMergeCut = groundMergeCut, ipad = ipad,
                      groundCutHeight = groundCutHeight, steepSlope = steepSlope, clothSize = clothSize,
                      doFilter = doFilter, filterRes = filterRes, filterN = filterN, 
                      clip.radius = clip.radius, clip.trajectory.distance = clip.trajectory.distance,
                      clip.x = clip.x, clip.y = clip.y,
                      draw.trajectory = draw.trajectory, 
                      selector = selector, dtm.path = dtm.path,
                      trafoMatrix.path = preTrafoMatrix.path,
                      rainFilter = rainFilter,
                      additionalSlices = additionalSlices,
                      exportSlice.upperLimit = exportSlice.upperLimit, 
                      exportSlice.lowerLimit = exportSlice.lowerLimit, #m
                      exportClippedLAS = exportClippedLAS,
                      dirPath = dirPath)
    
  }
  
  
  
  path_newSlice <- paste0(trafo.matchDirPath, "_total_ground_veg/", 
                          trafo.matchOldSet, trafo.matchSlice)
  
  
  if(trafo.matchOldSet != ""){
    cat("Calculating transformation matrix into old set\n")
    movingTrees1 <- Sys.time()
    
    
      
    
    cat(" - reading slice at old tree list from", basename(path_oldSlice), "\n")
    co <- capture.output(oldSlice <- readLAS(path_oldSlice))
    if(match.voxel.old > 0){
      oldSlice <- voxelize_points(oldSlice, match.voxel.old)
    }
    # keep the full old slice to match few points of new slice into it
    oldSlice@data$Z <- 0
    
    path_newSlice <- paste0(dirPath, "/_total_ground_veg/", fileFinder, trafo.matchSlice)
    if(!file.exists(path_newSlice)){
      stop(paste0("Could not transform, current slice file is missing in", 
                  path_newSlice, "\n"))
    }
    cat(" - reading current slice from", basename(path_newSlice), "\n")
    co <- capture.output(newSlice <- readLAS(path_newSlice))
    if(match.voxel.new > 0){
      newSlice <- voxelize_points(newSlice, match.voxel.new)
    }
    newSlice@data$Z <- 0
    
    old_mat <- as.matrix(oldSlice@data[, c("X", "Y", "Z")])
    old_mat <- old_mat[!duplicated.array(old_mat),]
    new_mat <- as.matrix(newSlice@data[, c("X", "Y", "Z")])
    new_mat <- new_mat[!duplicated.array(new_mat),]
    
    
    #plot3d(old_mat)
    
    cat(" - Merging lists via ICP - iterative closest point algorithm (Morpho)...\n")
    timeICP1 <- Sys.time()
    icp_result <- icpmat(new_mat, old_mat, mindist = 1, iterations = 100, type = "similarity")
    # moving = first, new_mat
    # fix = second, old_mat  
    timeICP2 <- Sys.time()
    cat("      Done in a ")
    print.difftime(round(timeICP2 - timeICP1,1))
    #cat("\n")
    #plot3d(old_mat, col = "blue")
    #plot3d(icp_result, col = "gold", add = T)
    
    #plot3d(old_mat, col = "blue")
    #plot3d(new_mat, col = "red", add = T)
    
    
    matr <- computeTransform(new_mat, icp_result)
    
    
    
    
    
    cat(" - final conversion matrix ")
    prmatrix(matr, rowlab=rep("  ",4), collab=rep("",4))
    
    
    groundPath <- v.env$groundPath
    matrixSavePath <- paste0(dirPath, groundPath, "trafo_from_", trafo.matchOldSet, "_to_",  fileFinder, ".txt")
    cat(" - saving matrix into", basename(matrixSavePath), "\n\n")
    write.table(matr, matrixSavePath, 
                row.names = F, sep = "\t", col.names = F)
    
    cat("Drawing final image to check for right transformation... ")
    
    drawSlice1 <- decimate_points(newSlice, random(300))
    drawSlice2 <- decimate_points(oldSlice, random(300))
    imagePath <- paste0(dirPath, "_images/transformed_vegetation/")
    if(!dir.exists(imagePath)) dir.create(imagePath, recursive = T)
    png(paste0(imagePath, fileFinder, "_into_", trafo.matchOldSet, ".png"), type = "cairo", 
        height = 4000, width = 4000)
    plot(drawSlice@data$X, drawSlice@data$Y, cex = 0.001, asp = 1, col = "black")
    points(drawSlice@data$X, drawSlice@data$Y, cex = 0.001, col = "red")
    dev.off()
    cat("done!\n")
    
    
    
    
    
    
    
  } 
  
  
  
  
  
  extractVegetation(LASfile = LASfile, fileFinder = fileFinder, 
                    groundMergeCut = groundMergeCut, ipad = ipad,
                    groundCutHeight = groundCutHeight, steepSlope = steepSlope, clothSize = clothSize,
                    doFilter = doFilter, filterRes = filterRes, filterN = filterN, 
                    clip.radius = clip.radius, clip.trajectory.distance = clip.trajectory.distance,
                    clip.x = clip.x, clip.y = clip.y,
                    draw.trajectory = draw.trajectory, 
                    selector = selector, dtm.path = dtm.path,
                    trafoMatrix.path = matrixSavePath,
                    rainFilter = rainFilter,
                    additionalSlices = additionalSlices,
                    exportSlice.upperLimit = exportSlice.upperLimit, 
                    exportSlice.lowerLimit = exportSlice.lowerLimit, #m
                    exportClippedLAS = exportClippedLAS,
                    dirPath = dirPath)
  
  
  
  
  
  movingTrees2 <- Sys.time()
  cat("Moving old tree list completed in a ")
  print.difftime(round(movingTrees2 - movingTrees1,1))
  cat("\n\n")
  
  
}





#' Ground classification of raw point cloud
#'
#' Reads in a .las or .laz file from disk and extracts the vegetation cloud
#' using the cloth simulation filter (scf) via lasground from the lidR package
#'
#' @param LASfile input path of full (or cropped) geoslam .las data
#' @param fileFinder user defined name to find this certain dataset in further processing
#' @param groundCutHeight in m, lower than this is all classified as ground (default: 0.5 m)
#' @param steepSlope in mountainous areas set it TRUE for better ground model
#' @param clothSize resolution of the cloth grid
#' @param filterRes noise filter resolution for adjacent 3x3x3 = 27 voxels, 2 cm (strong filter) takes longer time than 5 cm (medium filter)
#' @param filterN number of neighbour points, "strength" of filter, how many neighbours there need to be in all 27 neighbour voxels, (8 would be a medium, 27 a strong filter)

#' @param exportSlice.upperLimit in m, upper limit of single slice exported after ground classification for previewing stems and tree detection, set to 0 to disable slice export
#' @param exportSlice.lowerLimit in m, lower limit of single exported slice, default = 1 m

#' @param additionalSlices should additional slices be saved as .laz files
#' @param groundMergeCut in cm, how high should additional slice from ground be for a merging dtm (multiple tiles)

#' @param clip.trajectory.distance can be set to cut the input file according to the trajectory
#' @param clip.radius can be set to cut a circle from 0 0 with the given radius of the big input file
#' @param exportClippedLAS if set TRUE, then a laz file will be exported according to clip.radius
#' @export
extractVegetation <- function(LASfile, fileFinder, groundMergeCut = 0, ipad = FALSE,
                              groundCutHeight = 1.0, steepSlope = TRUE, clothSize = 0.10,
                              doFilter = FALSE, filterRes = 0.05, filterN = 27, 
                              clip.radius = 0, clip.trajectory.distance = 0,
                              clip.x = 0, clip.y = 0,
                              draw.trajectory = TRUE, 
                              selector = "xyzcit0RGB", dtm.path = "",
                              trafoMatrix.path = "",
                              rainFilter = 0,
                              additionalSlices = TRUE,
                              exportSlice.upperLimit = 3, exportSlice.lowerLimit = 1, #m
                              exportClippedLAS = FALSE,
                              dirPath = paste0(getwd(), "/")){

  
  if(length(LASfile) > 1){
    stop("More than one LASfile provided!\n")
  }
  if(length(fileFinder) > 1){
    stop("More than one fileFinder provided!\n")
  }
  if(length(trafoMatrix.path) > 1){
    stop("More than one transformation file provided!\n")
  }
  if(trafoMatrix.path != ""){
    if(!file.exists(trafoMatrix.path)){
      stop(paste0("Transformation file", trafoMatrix.path, "not found!\n"))
    }
  }
  if(dtm.path != ""){
    
    if(!file.exists(dtm.path)){
    stop(paste0("DTM file in", dtm.path, "not found!\n"))
      }
    }
  
  allStart <- Sys.time()
  fileFinder <- removeUmlaut(fileFinder)
  if(ipad){
    selector <- "xyzcitRGB0"
  }
  ##DEBUGGING
  if(2 == 1){
    groundCutHeight = 1
    steepSlope = TRUE
    clothSize = 0.10

    doFilter = FALSE
    filterRes = 0.05
    filterN = 27

    clip.radius = 0
    exportSlice.upperLimit = 3
    exportSlice.lowerLimit = 1 #m
    additionalSlices = T

    groundPath <- paste0("_total_ground_veg/")
  }


  groundPath <- v.env$groundPath
  if(!dir.exists(paste0(dirPath, groundPath))) dir.create(paste0(dirPath, groundPath))
  imgPath <- "_images/"
  if(!dir.exists(paste0(dirPath, imgPath))) dir.create(paste0(dirPath, imgPath))


  if(clip.trajectory.distance > 0 && exportClippedLAS){
    out.path.clip <- paste0(dirPath, "_traj", clip.trajectory.distance, "m/")
    if(!dir.exists(out.path.clip)){
      dir.create(out.path.clip)
      cat("Created the clipped directory", out.path.clip, "\n")
    }
  }

  if(clip.radius > 0 ){


    if(exportClippedLAS){
      out.path.clip <- paste0(dirPath, "_in_radius_", clip.radius, "m/")
      if(!dir.exists(out.path.clip)){
        dir.create(out.path.clip)
        cat("Created the clipped directory", out.path.clip, "\n")
      }
    }


     slicePath <- groundPath



    if(!dir.exists(slicePath)) dir.create(slicePath)

  }


  sink(paste0(dirPath, groundPath, fileFinder, "_extractVegetation_", format(Sys.time(), "%Y%m%d_%H%M"), "_Rcons.txt"), append = TRUE, split = TRUE)
  cat("Creating a ground and vegetation files for set", fileFinder, "\n")
  cat("Today is", format(Sys.time()), "\n")
  if(clip.trajectory.distance){
    cat(paste0("Clipping input file to trajectory +", clip.trajectory.distance, "m.\n"))
  }
  if(dtm.path != ""){
    cat("Normalizing with ground model from", basename(dtm.path),"\n")
  }
  if(trafoMatrix.path != "") cat("Using transformation file", paste0("\"", basename(trafoMatrix.path), "\""),"\n")
  cat("Working at", dirPath, "\n\n")
  #if(!dir.exists(dirPath)) dir.create(dirPath)
  
  

  allStart <- Sys.time()

  tempName <- basename(LASfile)
  
  
  
  
  # TRAJECTORY SEARCH
  tryCatch({
    
    finLaz <- strfind(tempName, "_100pct")
    if(is.null(finLaz)){
      finLaz <- strfind(tempName, ".las")
    }
    if(is.null(finLaz)){
      finLaz <- strfind(tempName, ".laz")
    }
    cat("Searching trajectory in",paste0(dirname(LASfile), "/", substr(tempName, 1, finLaz-1),"...\n    "))
    
    trajfile <- paste0(dirname(LASfile),
                       "/", substr(tempName, 1, finLaz-1), "_results_traj.txt")
    if(file.exists(trajfile)) cat("..._results_traj.txt ")
    if(!file.exists(trajfile)){
      trajfile <- paste0(dirname(LASfile),
                         "/", substr(tempName, 1, finLaz-1), ".gs-traj")
      if(file.exists(trajfile)) cat("...x.gs-traj ")
    }
    if(!file.exists(trajfile)){
      trajfile <- paste0(dirname(LASfile),
                         "/", substr(tempName, 1, finLaz-1), ".txt")
      if(file.exists(trajfile)) cat("...x.txt ")
    }
    if(!file.exists(trajfile)){
      trajfile <- paste0(dirname(LASfile),
                         "/", substr(tempName, 1, finLaz-1), "_traj.txt")
      if(file.exists(trajfile)) cat("..._traj.txt ")
    }
    if(!file.exists(trajfile) & clip.trajectory.distance > 0){
      cat("No trajectory found!\n")
      return("Clipping cannot be done without trajectory...\n\n")
    }
    
    # copy trajectory
    oneDone <- FALSE
    txtExists <- FALSE
    plyExists <- FALSE
    
    if(file.exists(trajfile)){
      txtExists <- TRUE
      file.copy(trajfile, paste0(dirPath, groundPath, fileFinder, "_traj.txt"), overwrite = T)
      cat("ok!    ")
    } else {
      cat("xXx non ")
    }
    oneDone <- TRUE
    cat("..._results_traj_time.ply ")
    
    trajfile <- paste0(dirname(LASfile), "/", substr(tempName, 1, finLaz-1), "_results_traj_time.ply")
    if(file.exists(trajfile)){
      plyExists <- TRUE
      file.copy(trajfile, paste0(dirPath, groundPath, fileFinder, "_traj.ply"), overwrite = T)
      cat("ok!\n")
    } else {
      trajfile <- paste0(dirname(LASfile), "/", substr(tempName, 1, finLaz-1), ".ply")
      if(file.exists(trajfile)){
        plyExists <- TRUE
        file.copy(trajfile, paste0(dirPath, groundPath, fileFinder, "_traj.ply"), overwrite = T)
        cat("ok!    ")
      } else {
        cat("xXx non ")
      }
    }
    
    cat("\n")
    if(txtExists && plyExists){
      cat("2 of 2 files sucessfully copied, find them as",
          paste0(fileFinder, "_traj.txt&ply\n"))
    } else{
      if(txtExists){
        cat("1 of 2 files sucessfully copied, find them as",
            paste0(fileFinder, "_traj.txt\n"))
      }
      if(plyExists){
        cat("1 of 2 files sucessfully copied, find them as",
            paste0(fileFinder, "_traj.ply\n"))
      }
    }
    cat("\n")
    
  }, error = function(error_condition) {
    cat("-> problem with the trajectory! ")
    if(oneDone) cat("Only the .txt file was copied!\n")
  })
  
  traj <- NA
  if(clip.trajectory.distance > 0 || draw.trajectory && txtExists){
    clipTime <- Sys.time()
    cat("Reading trajectory... ")
    try({
      traj <- read.csv(paste0(dirPath, groundPath, fileFinder, "_traj.txt"), sep = " ")
      if(is.element("X.time", colnames(traj))){
        duration_sec <- max(traj$X.time) - min(traj$X.time)
      } else if(is.element("X..world_time", colnames(traj))){
        duration_sec <- max(traj$X..world_time) - min(traj$X..world_time)
      } else {
        duration_sec <- max(traj[, 1]) - min(traj[, 1])
      }
      cat("(scan went", round(duration_sec/60,1), "mins)\n")
    })
    
    # GeoSLAM trajectory contains approx. 83 points per second,
    # we reduce them to 2 points per second for alpha-hulling
    #traj$X.time[2] - traj$X.time[1] original spacing of points is 1/100 of a second
    #traj2$X.time[2] - traj2$X.time[1] # now trajectory is spaced every 0.5 seconds
    
    traj2 <- traj[seq(from = 1, to = nrow(traj), by = 50),]
    traj <- traj2
    traj$col <- rainbow(length(traj[,1]), end = 0.7, rev = T)
    rm(traj2)
  }
  
  #if(tooBig){
  #  cat("tooBig: Reading only every 3rd point: ", LASfile, "...\n", sep = "")
  #  co <- capture.output(big <- readLAS(paste0(LASfile), select = selector, filter = "-keep_every_nth 3"))
  #  if(clothSize == 0.1){
  #    clothSize <- 0.5
  #  }
  # READING INPUT FILES SECTION ####
  {
    cat("Reading in: ", basename(LASfile), "... \n", sep = "")
    cat("(from ", dirname(LASfile), ")...", sep = "")
    readTime1 <- Sys.time()
    if(clip.radius > 0){
      cat("\nReading only a circle of r =",clip.radius,"m")
      if(clip.x != 0 || clip.y != 0){
        cat(paste0(" - base point is set to x=", round(clip.x, 2), " and y=", round(clip.y, 2), "m"))
      }
      cat("... ")
      co <- capture.output(big <- readLAS(LASfile, 
                                          select = selector, 
                                          filter = paste0("-keep_circle ", clip.x, " ", 
                                                          clip.y, " ",clip.radius)))
      
      
      oldHeader <- readLASheader(LASfile)
      pointsBefore <- oldHeader@PHB$`Number of point records`
      pointsAfter <- big@header@PHB$`Number of point records`
      pointsLost <- pointsBefore - pointsAfter
      
      cat("done!\n")
      cat("Remain", thMk(pointsAfter), "pts in the", paste0("r=", clip.radius, " m"), "circle.\n")
      cat("We lost", thMk(pointsLost), "pts (or", round(pointsLost/pointsBefore*100,1), "% of original pts)\n")
    } else {
      
      co <- capture.output(big <- readLAS(LASfile, select = selector))
      pointsBefore <- big@header@PHB$`Number of point records`
      pointsAfter <- pointsBefore
      pointsLost <- 0L
      cat("done!\n")
    }
    readTime2 <- Sys.time()
    cat("Reading input files took a ")
    print.difftime(round(readTime2 - readTime1, 1))
    cat("\n")
  }


  if(!is.element("Intensity", colnames(big@data))){
    big@data$Intensity <- 5L
    warning("NO INTENSITY VALUES PRESENT!!!\n\n")
  }
  cat("Checking if intensity values are present...\n")
  if(max(big@data$Intensity) == 0){
    big@data$Intensity <- 65535L
#    writeLAS(big, "D:/nowlas.laz")
  }
  if(max(big@data$Intensity) <= 256L){
    cat("Expanding intensity from 255 to 65535...\n")
    big@data$Intensity <- as.integer(round((big@data$Intensity / 255L) * 65535))
    #    writeLAS(big, "D:/nowlas.laz")
  }


  if(big@header@PHB$`X scale factor` < 0.01) big@header@PHB$`X scale factor` <- 0.00001 # correcting scale factors
  if(big@header@PHB$`Y scale factor` < 0.01) big@header@PHB$`Y scale factor` <- 0.00001
  if(big@header@PHB$`Z scale factor` < 0.01) big@header@PHB$`Z scale factor` <- 0.00001
  if(big@header@PHB$`X offset` < 0) big@header@PHB$`X offset` <- 0 # correcting offsets, else cannot write las with negative offset... stupid rules...
  if(big@header@PHB$`Y offset` < 0) big@header@PHB$`Y offset` <- 0
  if(big@header@PHB$`Z offset` < 0) big@header@PHB$`Z offset` <- 0
  if(clip.radius > 0 && exportClippedLAS){
    writeLAS(big, paste0(out.path.clip, fileFinder, "_", clip.radius, "m.laz"))
  }

  cat("There were", thMk(big@header@PHB$`Number of point records`), "points found.\n")
  #lidR::plot(big)

  if(trafoMatrix.path != "" & file.exists(trafoMatrix.path)){
    library("Morpho")
    cat("\nTransforming point cloud with file", basename(trafoMatrix.path),"\n")
    tryCatch({
      tf1 <- Sys.time()
      suppressWarnings(trafoMat <- as.matrix(read.table(trafoMatrix.path)))
      trafoMat <- trafoMat[c(1:4), c(1:4)]
      prmatrix(trafoMat, rowlab=rep("",4), collab=rep("",4))
      cat("   -> convert, ")
      tempCoords <- applyTransform(as.matrix(big@data[, c(1:3)]), trafo = trafoMat)
      cat("assign, ")
      big@data$X <- round(tempCoords[,1], 6)
      big@data$Y <- round(tempCoords[,2], 6)
      big@data$Z <- round(tempCoords[,3], 6)
      cat("clean, ")
      rm(tempCoords)
      gc()

      cat("done in a ")
      print.difftime(round(Sys.time() - tf1, 1))
      cat("\n")
    }, error = function(error_condition) {
      cat("NO TRANSFORMATION WAS DONE!!!\n\n\n")
      warning("Transformation process was not sucessful!")
    })
  }



  if(trafoMatrix.path != "" & file.exists(trafoMatrix.path) & txtExists){
    cat("\n\nTransforming trajectory with file", basename(trafoMatrix.path),"\n")
    tryCatch({
      ttf1 <- Sys.time()
      cat("Reading trajectory for transformation... \n")
      traj <- read.csv(paste0(dirPath, groundPath, fileFinder, "_traj.txt"), sep = " ")
      suppressWarnings(trafoMat <- as.matrix(read.table(trafoMatrix.path)))
      trafoMat <- trafoMat[c(1:4), c(1:4)]
      cat("convert, ")
      tempCoords <- applyTransform(x = as.matrix(traj[, c(2:4)]), trafo = trafoMat)
      cat("assign, ")
      traj$x <- round(tempCoords[,1], 6)
      traj$y <- round(tempCoords[,2], 6)
      traj$z <- round(tempCoords[,3], 6)
      cat("save, ")
      colnames(traj)[1] <- "%time"
      write.table(traj, paste0(dirPath, groundPath, fileFinder, "_traj.txt"),
                  sep = " ", row.names = F, quote = F, na = "")
      cat("clean, ")
      rm(tempCoords)
      gc()

      cat("done in a ")
      print.difftime(round(Sys.time() - ttf1, 1))
      
      traj2 <- traj[seq(from = 1, to = nrow(traj), by = 50),]
      traj <- traj2
      traj$col <- rainbow(length(traj[,1]), end = 0.7, rev = T)
      rm(traj2)
      
      cat("")
    }, error = function(error_condition) {
      cat("NO TRANSFORMATION WAS DONE!!!\n\n\n")
      warning("Transformation process was not sucessful!")
    })
  }


  cat("\n")


  
  ## TRAJECTORY PNG IF CIRCLE OR FULL AREA ####
  if(clip.trajectory.distance == 0){
    if(clip.radius > 0){
      big_sm <- decimate_points(big, random(10))
      tempHeight <- clip.radius*20
      if(tempHeight < 500) tempHeight <- 500
      png(paste0(dirPath, imgPath, fileFinder, "_circle_traj.png"),
          height = tempHeight, width = tempHeight, type = "cairo")
      yLims <- c(clip.y - clip.radius, clip.y + clip.radius)
      xLims <- c(clip.x - clip.radius, clip.x + clip.radius)
      nowTitle <-  paste0(fileFinder, " circle radius ", clip.radius, 
                          " m (area=", round(clip.radius^2*pi/10000,3), "ha)")
    } else {
      big_sm <- decimate_points(big, random(1))
      
      tempHeight <- (big_sm@header@PHB$`Max Y` - big_sm@header@PHB$`Min Y`)*5
      tempWidth <-  (big_sm@header@PHB$`Max X` - big_sm@header@PHB$`Min X`)*5
      if(tempHeight < 500 || tempWidth < 500){
        if(tempHeight < tempWidth){
          tempWidth <- tempWidth/tempHeight * 500
          tempHeight <- 500
        } else {
          tempHeight <- tempHeight/tempWidth * 500
          tempWidth <- 500
        }
      } 
      png(paste0(dirPath, imgPath, fileFinder, "_plot_traj.png"),
          height = tempHeight, width = tempWidth, type = "cairo")
      yLims <- c(big_sm@header@PHB$`Min Y`, big_sm@header@PHB$`Max Y`)
      xLims <- c(big_sm@header@PHB$`Min X`, big_sm@header@PHB$`Max X`)
      nowTitle <-  paste0(fileFinder, " full plot")
    }
    plot(big_sm$Y ~ big_sm$X, cex = 0.0001, asp = 1,
         xlab = "x [m]", ylab = "y [m]", xlim = xLims, ylim = yLims,
         main = nowTitle)
    if(clip.radius > 20){
      draw.circle(clip.x, clip.y, 20, lwd = 1, lty = 2, border = "red")
    } else if(clip.radius > 10){
      draw.circle(clip.x, clip.y, 10, lwd = 1, lty = 2, border = "red")
    }
    if(length(traj) > 1){
      points(traj$y ~ traj$x, col = traj$col, cex = 0.5, pch = 16)
      legend("topright", col = c("blue", "red"), lwd = 3,
             legend = c("start", paste0(round(duration_sec/60,1), " min")),
             bty = "n")
    }
    #legend("bottomleft", legend = c("1", "2", "3"), cex = 0.5)
    legend("bottomleft", legend = c(paste0(" ", thMk(pointsBefore), " pts"),
                                    paste0("-",thMk(pointsAfter), " pts"),
                                    paste0("", strcat(rep("-", nchar(thMk(pointsAfter))*2)), "---"),
                                    paste0("  ", thMk(pointsLost), " pts lost")), cex = 1.1, bty = "n")
    
    rm(big_sm)
    #plot(hull.traj, add=TRUE, border= 2, lwd = 2)
    dev.off()
  }
  
  
  
      
      
      
      
      
      


  # TRAJECTORY CLIPPING SECTION ####

  if(clip.trajectory.distance > 0 && txtExists){
    cat("Clipping total cloud to trajectory +", clip.trajectory.distance, "m radius (a=25m)...")
    
    # shifting trajectory for the sake of polygon hulling, goes wrong with gps coordinates!
    #dups <- duplicated.data.frame(data.frame("x" = traj$x, "y" = traj$y))
    shiftX <- min(traj$x)
    shiftY <- min(traj$y)
    if(shiftX > 1000){
      cat(" shft x", shiftX, "m")
      traj$x <- traj$x - shiftX
    }
    if(shiftY > 1000){
      cat(" shft y", shiftY, "m")
      traj$y <- traj$y - shiftY
    }
    
    borderHull <- ashape(traj$x, traj$y, alpha = 25)
    owin.window <-owin(xrange=range(traj$x), yrange=range(traj$y))
    trees_edges_ppp <- psp(borderHull$edges[,3], borderHull$edges[,4], borderHull$edges[,5], borderHull$edges[,6], window=owin.window)
    hull.traj <- dilation(trees_edges_ppp, r=clip.trajectory.distance)
    owin.traj <- owin(poly=hull.traj$bdry[[1]])
    area.traj <- area.owin(owin.traj)

    big_sm <- decimate_points(big, random(1))

    outer <- Polygon(hull.traj$bdry[[1]]) # for clipping LAS point cloud

    if(shiftX > 1000){
      traj$x <- traj$x + shiftX
      outer@labpt <- outer@labpt + c(shiftX, 0)
      outer@coords[,1] <- outer@coords[,1] + shiftX
      hull.traj$xrange <- hull.traj$xrange + shiftX
    }
    if(shiftY > 1000){
      traj$y <- traj$y + shiftY
      outer@labpt <- outer@labpt + c(0, shiftY)
      outer@coords[,2] <- outer@coords[,2] + shiftY
      hull.traj$yrange <- hull.traj$yrange + shiftY
    }


    pointsBefore <- big@header@PHB$`Number of point records`
    big <- clip_roi(big, outer)
    pointsAfter <- big@header@PHB$`Number of point records`

    pointsLost <- pointsBefore - pointsAfter
    cat("done!\n")
    cat("Remain", thMk(pointsAfter), "pts after clipping to shape.\n")
    cat("We lost", thMk(pointsLost), "pts (or", round(pointsLost/pointsBefore*100,1), "% of original pts)\n")

    ### TRAJECTORY CLIP PNG ####
    png(paste0(dirPath, imgPath, fileFinder, "_traj_clipping.png"),
        height = diff(hull.traj$yrange)*5, width = diff(hull.traj$xrange)*5, type = "cairo")
    plot(big_sm$Y ~ big_sm$X, cex = 0.0001, asp = 1,
         xlab = "x [m]", ylab = "y [m]", xlim = hull.traj$xrange, ylim = hull.traj$yrange,
         main = paste0(fileFinder, "_traj.txt (+", clip.trajectory.distance, "m) area=", round(area.traj/10000,2), "ha"))
    points(traj$y ~ traj$x, col = traj$col, cex = 0.5, pch = 16)
    legend("topright", col = c("blue", "red"), lwd = 3,
           legend = c("start", paste0(round(duration_sec/60,1), " min")),
           bty = "n")
    #legend("bottomleft", legend = c("1", "2", "3"), cex = 0.5)
    legend("bottomleft", legend = c(paste0(" ", thMk(pointsBefore), " pts"),
                                    paste0("-",thMk(pointsAfter), " pts"),
                                    paste0("", strcat(rep("-", nchar(thMk(pointsAfter))*2)), "---"),
                                    paste0(" ", thMk(pointsLost), " pts lost")), cex = 1.1, bty = "n")

    lines(outer@coords, col = 2, lwd = 2)
    #plot(hull.traj, add=TRUE, border= 2, lwd = 2)
    dev.off()
    clipTime2 <- Sys.time()
    cat("Trajectory clipping took a ")
    print.difftime(round(clipTime2 - clipTime,1))


    rm(big_sm)
    if(exportClippedLAS){
      writeLAS(big, paste0(out.path.clip, fileFinder, "_traj", clip.trajectory.distance, "m.laz"))
    }
    gc()
  }

  
  # DTM SECTION ####

  ## DTM GIVEN ####
  if(file.exists(dtm.path)){
    # if DTM is given, we only need to cut it once to absolute values (need no puffer)

    cat("\nDTM-Section:\n")
    gstart <- Sys.time()
    tryCatch(
      {
        # read in raster file
        cat("Reading input dtm from",dtm.path,"\n")
        dtm_y <- raster(dtm.path)

        dtSave <- crop(dtm_y, extent(c(big@header$`Min X`-0.5, big@header$`Max X`+0.5,
                                       big@header$`Min Y`-0.5, big@header$`Max Y`+0.5)))
        writeRaster(dtSave, paste0(dirPath, groundPath, fileFinder, "_ground_clip.tif"), overwrite = TRUE)
        writeRaster(dtSave, paste0(dirPath, groundPath, fileFinder, "_ground_clip.grd"), overwrite = TRUE)
        rm(dtSave)
      }, error = function(error_condition) {
        warning("Error in reading the dtm-model, file", basename(dtm.path), "not found!")
        return()
      })



    ## DEM - same res as DTM and 0.01 disk size ####
    dem.res <- raster::res(dtm_y)[1]
    dem.p2r <- 0.01

    cat(paste0("(*) CANOPY DEM res",dem.res*100,"cm point radius",dem.p2r*100,"cm... "))
    co <- capture.output(canop <- grid_canopy(big, res = dem.res, p2r(dem.p2r)))
    #plot(canop, main = paste0(groundCutHeight," DEM res",dem.res, " p2r",dem.p2r))

    if(rainFilter > 0){
      rf1 <- Sys.time()
      for(i in 1:rainFilter){
        if(i == 1){
          cat("RAIN FILTER ")
          normi <- big
        }
        cat(i)
        cat("x ")
        normi <- normalize_height(normi, canop)
        normi <- filter_poi(normi, Z <= -0.01)
        normi <- unnormalize_height(normi)
        co <- capture.output(canop <- grid_canopy(normi, res = dem.res, p2r(dem.p2r)))
        #plot(canop, main = paste0(groundCutHeight," DEM res",dem.res, " p2r",dem.p2r))
      }
      rf2 <- Sys.time()
      timeRainFilter <- as.difftime(rf2 - rf1)
      cat("DONE IN", round(timeRainFilter,1), units(timeRainFilter),"")
    }
    writeRaster(canop, paste0(dirPath, groundPath, fileFinder, "_dem.tif"), overwrite = TRUE)
    writeRaster(canop, paste0(dirPath, groundPath, fileFinder, "_dem.grd"), overwrite = TRUE)
    cat(" saved!\n")


    # NDOM ####
    tryCatch(
      {
        #dtCrop <- crop(dtm_y, extent(c(big@header$`Min X`, big@header$`Max X`, big@header$`Min Y`, big@header$`Max Y`)))
        dtCrop <- dtm_y
        cat(paste0("(*) NDOM in assembly... "))
        ndom <- canop - dtCrop
        ndom@extent
        #ndom@data@values[is.na(ndom@data@values)] <- 0
        ndom@data@values[!is.na(ndom@data@values) & ndom@data@values < 0] <- 0
        ndom@data@values[!is.na(ndom@data@values) & ndom@data@values > 44] <- NA
        ndom@data@values[!is.na(ndom@data@values) & ndom@data@values > 41] <- 41
        if(clip.radius != 0){
          #create a 40 m radius clipping
          cc_20 <- dismo::circles(data.frame(clip.x, clip.y), lonlat = F, clip.radius)
          ndom <- mask(crop(ndom, polygons(cc_20)),polygons(cc_20))
        }
        ndom@data@values[!is.na(ndom@data@values)][1:100] <- 41
        writeRaster(ndom, paste0(dirPath, groundPath, fileFinder, "_ndom.tif"), overwrite = TRUE)
        writeRaster(ndom, paste0(dirPath, groundPath, fileFinder, "_ndom.grd"), overwrite = TRUE)
        cat(" saved!\n")
        ndom@data@values[is.na(ndom@data@values)] <- 0

        png(filename = paste0(dirPath, imgPath, fileFinder,"_NDOM_raw.png"), 
            width = 800, height = 800, type = "cairo")
        par(xpd = F, mar = c(2,0,5,0), oma = c(0,0,0,0), xaxs='i', yaxs='i')
        if(clip.radius != 0){
          plot(0, type = "n",
               xlim = c(-clip.radius + clip.x, clip.radius + clip.x),
               ylim = c(-clip.radius + clip.y, clip.radius + clip.y), asp = 1, main = fileFinder, cex.main = 4, axes = F)
          raster::plot(ndom,  col = c(rainbow(19, start = 0.76, end = 0.60, rev = T), "#56008F", "#000000"), # every color is 1 m step (40 m span, 41 = taller than 40)
                       box = F, axes = F, legend = F,
                       asp = 1, add = T)
        } else {
          raster::plot(ndom,  col = c(rainbow(19, start = 0.76, end = 0.60, rev = T), "#56008F", "#000000"), # every color is 1 m step (40 m span, 41 = taller than 40)
                       box = F, axes = F, legend = F, cex.main = 4,
                       asp = 1, main = fileFinder)
        }
        axis(1, cex.axis = 2, lwd = 2)

        lgd_ <- c(">",seq(40, 1, by = -2))
        lgd_[nchar(lgd_) < 2] <- paste0("  ", lgd_[nchar(lgd_) < 2])
        legend("topright", title = "NDOM [m]",
               legend = lgd_, pt.cex = 3,
               fill = c("#000000", "#56008F", rainbow(19, start = 0.76, end = 0.60, rev = F)),
               border = NA,
               y.intersp = 0.7, cex = 1, text.font = 1, box.lty = 0)
        if(clip.radius != 0){
          clip(x1 = -clip.radius + clip.x + 0.5, x2 = clip.radius + clip.x - 0.5,
               y1 = -clip.radius + clip.y + 0.5, y2 = clip.radius + clip.y - 0.5)
        }
        abline(h = clip.y, v = clip.x, lwd = 1)
        if(clip.radius > 20){
          draw.circle(clip.x, clip.y, 20, lwd = 2, lty = 2)
        } else if(clip.radius > 10){
          draw.circle(clip.x, clip.y, 10, lwd = 2, lty = 2)
        }
        dev.off()

      }, error = function(error_condition) {
        warning("NDOM not created... Is that a problem?\n")
      })


    gstop <- Sys.time()
    cat("ONLY DEM and NDOM done for provided DTM in a")
    print.difftime(round(gstop - gstart, 1))



    # normalization and slice cutting ####
    tryCatch(
      {
        cat("Normalizing height with input model... ")
        big <- normalize_height(big, dtCrop, na.rm = TRUE) # need to save it in that intermediate object or it cannot unnormalize anymore
        cat("done!\n")
      }, error = function(error_condition) {
        cat("Error in using the dtm-model...\n")
      })

    cat("\nGround cutoff at height", groundCutHeight * 100, "cm\n")
    ground <- filter_poi(big, Z <= groundCutHeight)
    ground@data$Classification <- 2L
    ground <- unnormalize_height(ground)
    groundPts <- ground@header@PHB$`Number of point records`

    vegetation <- filter_poi(big, Z > groundCutHeight)
    vegetation <- unnormalize_height(vegetation)
    vegetation@data$Classification <- 1L

    vox.gr <- 0.05
    cat("Voxelize all ground points with a raster of",vox.gr*100,"cm... ")
    {
      t2 <- Sys.time()
      gr_vox <- TreeLS::tlsSample(ground, TreeLS::smp.voxelize(vox.gr))
      t3 <- Sys.time()
      cat("done by treeLS. ")
      groundFullVoxPts <- gr_vox@header@PHB$`Number of point records`
      writeLAS(gr_vox, paste0(dirPath, groundPath, fileFinder, "_ground_",groundCutHeight*100,"cm_vox.laz"))
      print.difftime(round(t3-t2,1))
    }

    # writing laz files, sor filtering ####
    {
      cat("Writing out ground file ", paste0(dirPath, groundPath, fileFinder, "_ground.laz... \n"))
      writeLAS(ground, file = paste0(dirPath, groundPath, fileFinder, "_ground.laz"))
      LAS_ground <<- ground
      # groundSm <- decimate_points(ground, random(200))
      # groundSmPts <- groundSm@header@PHB$`Number of point records`
      # cat(" and", paste0(fileFinder, "_ground_sm.las... "))
      # writeLAS(groundSm, file = paste0(dirPath, groundPath, fileFinder, "_ground_sm.las"))
      #
      # rm(groundSm)
      if(!additionalSlices) rm(ground) # keep the ground for additional slices
      gc()

      # "SOR" Noise filtering
      if(doFilter){
        cat("\nStarting with noise filtering (classify_noise with ivf algorithm in lidR):\n")
        cat("Voxel resolution:", filterRes, "\nminimum number of points in these voxels:", filterN, "\n")
        startSOR <- Sys.time()
        vegetation <- classify_noise(las = vegetation, algorithm = ivf(res = filterRes, n = filterN))
        stopSOR <- Sys.time()

        noiseP <- (sum(vegetation@data$Classification == 18))
        allP <- vegetation@header@PHB$`Number of point records`
        cat("Classified", thMk(noiseP), "points as noise, that is", round(noiseP/allP*100, 2), "%.\n")

        cat("Noise filtering done in a ")
        print.difftime(round(stopSOR - startSOR, 1))
        cat("\nWriting out")
      } else {
        cat("    and")
      }

      cat(" vegetation file", paste0(dirPath, groundPath, fileFinder, "_raw_veg.laz... "))
      writeLAS(vegetation, file = paste0(dirPath, groundPath, fileFinder, "_raw_veg.laz"))
      LAS_veg <<- vegetation
      LAS_veg_name <<- fileFinder
      cat("done!\n")

      {
        cat("\n The ground file up to",groundCutHeight*100,"cm contains", thMk(groundPts), "points, \n")
        cat("       the 5 cm voxelized file contains",thMk(groundFullVoxPts),"points, \n")
        #cat("     and minimal voxel file up to 10 cm",thMk(groundMinVoxPts),"points.\n")
        cat("The raw vegetation file above contains", thMk(vegetation@header@PHB$`Number of point records`), "points.\n")
        cat("COMPARE: Comlete input cloud contains", thMk(big@header@PHB$`Number of point records`), "points (")
        cat(round(groundPts/big@header@PHB$`Number of point records`*100,0),"% ground)\n")
        if(1 == 2){
          #lidR::plot(vegetation)

          ints <- vegetation@data$Intensity
          cat("Intensity is scoping from",thMk(min(ints)),"to",thMk(max(ints)),"\n")
          #manual Histogram
          #hist(ints, freq = FALSE, breaks = c(0,(1:5)*10000,max(ints)))
          hist(ints, freq = FALSE)
          #plot(sort(ints[1:1000], decreasing = TRUE))
        }
      }

      rm(vegetation, ground, gr_vox)
      gc()

    }



    # create practical slices for easier reading data ####
    if(exportSlice.upperLimit > 0){
      allStop <- Sys.time()
      cat("\nGround and vegetation splitting completed in a ")
      print.difftime(round(allStop - allStart,1))

      cat("\nCreating additional slices: \n")



      slicePath <- paste0(dirPath, groundPath, fileFinder,
                          "_clusterSlice_100to300.laz")
      #voxSlicePath <- paste0(dirPath, groundPath, fileFinder,
      #                       "_clusterSlice_100to300_vox.laz")
      slice <- filter_poi(big, Z < 3,
                          Z > 1)
      slice_un <- unnormalize_height(slice)


      cat("   o) slice for tree detection from 100 to 300 cm... ")
      writeLAS(slice_un, slicePath)
      cat("done!\n")


      if(exportSlice.lowerLimit != 1 | exportSlice.upperLimit != 3){

        slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",
                            exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,".laz")
        #voxSlicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",
        #                       exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,"_vox.laz")


        #rm(dtm, dtm_fine, dtm3)
        slice <- filter_poi(big, Z < exportSlice.upperLimit, Z > exportSlice.lowerLimit)
        slice_un <- unnormalize_height(slice)

        cat("   o) desired special slice ", exportSlice.lowerLimit*100, "-", exportSlice.upperLimit*100, " cm... ", sep = "")
        writeLAS(slice_un, slicePath)
        cat("done!\n")
      }







      if(additionalSlices){
        cat("   o) additional slices... 120-140... ")
        # 0 to 2
        slice <- filter_poi(big, Z < 1.4, Z > 1.2)
        slice_un <- unnormalize_height(slice)
        slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_120to140.laz")
        writeLAS(slice_un, slicePath)

        cat("50-200... ")
        # 0 to 2
        slice <- filter_poi(big, Z < 2, Z > 0.5)
        slice_un <- unnormalize_height(slice)
        slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_50to200.laz")
        writeLAS(slice_un, slicePath)

        cat("300-500... ")
        slice <- filter_poi(big, Z < 5, Z > 3)
        slice_un <- unnormalize_height(slice)
        slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_300to500.laz")
        writeLAS(slice_un, slicePath)


        if(groundMergeCut > 0){
          cutter <- (groundMergeCut/100)
          cat(paste0("0-", groundMergeCut, "... "))
          slice <- filter_poi(big, Z < cutter)
          slice_un <- unnormalize_height(slice)
          slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_0_ground_", groundMergeCut, ".laz")
          writeLAS(slice_un, slicePath)
        }


        cat("done!\n")

        slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,".laz")
      }
      rm(slice, slice_un)
    }
  }  else {
    ## old stuff - creating new ground model
    # CREATING NEW DTM MODELS ####
    {

      cat("\nCreating new DTM:\n")
      gstart <- Sys.time()

      steepCliffsOfStone <- FALSE
      if(steepCliffsOfStone){
        fine_clothSize <- 0.02
        cat("Steep stone cliffs very fine cloth resolution", fine_clothSize*100, "cm, steep slopes:", TRUE, ", rigidness = 3 (less rugged terrain)... takes a while...\n")
        tground <- Sys.time()
        big <- classify_ground(big, csf(class_threshold = 3,
                                        sloop_smooth = T, rigidness = 1,
                                        cloth_resolution = 0.1), last_returns = FALSE)
        ground <- filter_poi(big, Classification == 2)
        ttground <- Sys.time()
        cat("Required a ")
        print.difftime(round(ttground - tground, 1))

        cat("Writing out ground file ", paste0(dirPath, groundPath, fileFinder, "_ground_2cm_cloth.laz... \n"))
        writeLAS(ground, file = paste0(dirPath, groundPath, fileFinder, "_ground_2cm_cloth.laz"))


        vox.gr <- 0.1
        cat("Voxelize minimal ground points with a raster of",vox.gr*100,"cm... ")
        {
          t2 <- Sys.time()
          gr_vox <- TreeLS::tlsSample(ground, TreeLS::smp.voxelize(vox.gr))
          t3 <- Sys.time()
          cat("done by treeLS ")
          groundMinVoxPts <- gr_vox@header@PHB$`Number of point records`
          writeLAS(gr_vox, paste0(dirPath, groundPath, fileFinder, "_ground_min_vox.laz"))
          cat("in a ")
          print.difftime(round(t3-t2,1))
        }
        plot(gr_vox)




        # also not all that perfect
        res.dtm =2
        k.dtm = 100000
        p.dtm = 5
        cat(paste0("(*) MINIMAL DTM gcut10cm res",res.dtm," k",k.dtm," p",p.dtm," ... "))
        dtmTime <- Sys.time()
        co <- capture.output(dtm_z <- grid_terrain(gr_vox, res = res.dtm,
                                                   algorithm = knnidw(k = k.dtm,  p = p.dtm)))


        print.difftime(round(Sys.time() - dtmTime, 1))
        #plot(dtm_z, main = paste0("10cm gcut ",clothSize*100,"cm cloth, res",res.dtm, " k",k.dtm," p",p.dtm))
        writeRaster(dtm_z, paste0(dirPath, groundPath, fileFinder, "_ground_min_1m_special.tif"), overwrite = TRUE)
        writeRaster(dtm_z, paste0(dirPath, groundPath, fileFinder, "_ground_min_1m_special.grd"), overwrite = TRUE)

        # changing the RESOLUTION of the DTM Raster (resample) finer density (10 cm)
        dtm_y <- raster(nrow = dtm_z@nrows * 10, ncol = dtm_z@ncols * 10)
        dtm_y@extent <- dtm_z@extent
        dtm_y <- resample(dtm_z, dtm_y, method='bilinear')
        #plot(dtm_fine)
        writeRaster(dtm_y, paste0(dirPath, groundPath, fileFinder, "_ground_min_special.grd"), overwrite = TRUE)
        writeRaster(dtm_y, paste0(dirPath, groundPath, fileFinder, "_ground_min_special.tif"), overwrite = TRUE)
        plot(dtm_y)




      }


      ## MIN - only keep 10 cm high points ####
      cat("Ground cutoff at height 10 cm for minimal ground points.\n")
      fine_clothSize <- 0.05
      if(clip.radius == 0){
        cat(" - - Very big file - better set the fine cloth also to a 25 cm rough one...\n")
        #cat(" - - We try still with a fine cloth of 5 cm - change it in f1 line 144 code R!\n")
        fine_clothSize <- 0.25
      }
      cat("Finest cloth resolution", fine_clothSize*100, "cm, steep slopes:", FALSE, ", rigidness = 3 (less rugged terrain)... takes a while...\n")
      #big <- lasground(big, csf(class_threshold = 0.05, cloth_resolution = 0.2, rigidness = 1))
      tground <- Sys.time()
      big <- classify_ground(big, csf(class_threshold = 0.1,
                                      sloop_smooth = FALSE, rigidness = 3,
                                      cloth_resolution = fine_clothSize), last_returns = FALSE)
      ground <- filter_poi(big, Classification == 2)
      ttground <- Sys.time()
      cat("Required a ")
      print.difftime(round(ttground - tground, 1))

      vox.gr <- 0.05
      cat("Voxelize minimal ground points with a raster of",vox.gr*100,"cm... ")
      {
        t2 <- Sys.time()
        gr_vox <- TreeLS::tlsSample(ground, TreeLS::smp.voxelize(vox.gr))
        t3 <- Sys.time()
        cat("done by treeLS ")
        groundMinVoxPts <- gr_vox@header@PHB$`Number of point records`
        writeLAS(gr_vox, paste0(dirPath, groundPath, fileFinder, "_ground_min_vox.laz"))
        cat("in a ")
        print.difftime(round(t3-t2,1))
      }

      res.dtm = 1
      k.dtm = 1000
      p.dtm = 2
      cat(paste0("(*) MINIMAL DTM gcut10cm res",res.dtm," k",k.dtm," p",p.dtm," ... "))
      dtmTime <- Sys.time()
      co <- capture.output(dtm_z <- grid_terrain(gr_vox, res = res.dtm,
                                                 algorithm = knnidw(k = k.dtm,  p = p.dtm)))
      print.difftime(round(Sys.time() - dtmTime, 1))
      #plot(dtm_z, main = paste0("10cm gcut ",clothSize*100,"cm cloth, res",res.dtm, " k",k.dtm," p",p.dtm))
      writeRaster(dtm_z, paste0(dirPath, groundPath, fileFinder, "_ground_min_1m.tif"), overwrite = TRUE)
      writeRaster(dtm_z, paste0(dirPath, groundPath, fileFinder, "_ground_min_1m.grd"), overwrite = TRUE)

      # changing the RESOLUTION of the DTM Raster (resample) finer density (10 cm)
      dtm_y <- raster(nrow = dtm_z@nrows * 10, ncol = dtm_z@ncols * 10)
      dtm_y@extent <- dtm_z@extent
      dtm_y <- resample(dtm_z, dtm_y, method='bilinear')
      #plot(dtm_fine)
      writeRaster(dtm_y, paste0(dirPath, groundPath, fileFinder, "_ground_min.grd"), overwrite = TRUE)
      writeRaster(dtm_y, paste0(dirPath, groundPath, fileFinder, "_ground_min.tif"), overwrite = TRUE)


      cat("\nGround cutoff at height", groundCutHeight * 100, "cm\n")
      cat("Cloth resolution", clothSize*100, "cm, steep slopes:", steepSlope, ", rigidness = 3 (less rugged terrain)\n")
      big <- classify_ground(big, csf(class_threshold = groundCutHeight,
                                      sloop_smooth = steepSlope, rigidness = 3,
                                      cloth_resolution = clothSize), last_returns = FALSE)
      ground <- filter_poi(big, Classification == 2)
      groundPts <- ground@header@PHB$`Number of point records`
      cat("Voxelize all ground points with a raster of",vox.gr*100,"cm... ")
      {
        t2 <- Sys.time()
        gr_vox <- TreeLS::tlsSample(ground, TreeLS::smp.voxelize(vox.gr))
        t3 <- Sys.time()
        cat("done by treeLS. ")
        groundFullVoxPts <- gr_vox@header@PHB$`Number of point records`
        writeLAS(gr_vox, paste0(dirPath, groundPath, fileFinder, "_ground_",groundCutHeight*100,"cm_vox.laz"))
        print.difftime(round(t3-t2,1))
      }




      # ROUGH - 1 m, 20 k and 2 p ####
      res.dtm = 1
      k.dtm = 20
      p.dtm = 2
      cat(paste0("(*) ROUGH DTM gcut",groundCutHeight*100,"cm res",res.dtm," k",k.dtm," p",p.dtm," ... "))
      dtmTime <- Sys.time()
      co <- capture.output(dtm_rough <- grid_terrain(gr_vox, res = res.dtm,
                                                     algorithm = knnidw(k = k.dtm,  p = p.dtm)))
      cat("done in a ")
      print.difftime(round(Sys.time() - dtmTime, 1))
      #plot(dtm_rough, main = paste0(groundCutHeight," gcut ",clothSize*100,"cm cloth, res",res.dtm, " k",k.dtm," p",p.dtm))
      writeRaster(dtm_rough, paste0(dirPath, groundPath, fileFinder, "_ground_rough.tif"), overwrite = TRUE)
      writeRaster(dtm_rough, paste0(dirPath, groundPath, fileFinder, "_ground_rough.grd"), overwrite = TRUE)






      # FINE - 0.2 m, 500 k and 0.5 p ####
      res.dtm = 0.2
      k.dtm = 500
      p.dtm = 0.5
      cat(paste0("(*) FINE DTM gcut",groundCutHeight*100,"cm res",res.dtm," k",k.dtm," p",p.dtm," ... "))
      dtmTime <- Sys.time()
      co <- capture.output(dtm_fine <- grid_terrain(gr_vox, res = res.dtm,
                                                    algorithm = knnidw(k = k.dtm,  p = p.dtm)))
      cat("done in a ")
      print.difftime(round(Sys.time() - dtmTime, 1))
      #plot(dtm_fine, main = paste0(groundCutHeight," gcut ",clothSize*100,"cm cloth, res",res.dtm, " k",k.dtm," p",p.dtm))
      writeRaster(dtm_fine, paste0(dirPath, groundPath, fileFinder, "_ground_fine.tif"), overwrite = TRUE)
      writeRaster(dtm_fine, paste0(dirPath, groundPath, fileFinder, "_ground_fine.grd"), overwrite = TRUE)


      if(1==2){
      #if(denoiseDEM){
        # would be nice, but is way to slow (takes up to 10 hours for 0.5 ha)
        big <- classify_noise(big, sor(k = 6, m = 2))
        noiseLAS <- filter_poi(big, Classification == LASNOISE)
        writeLAS(noiseLAS, file = paste0(dirPath, groundPath, fileFinder, "_noise.laz"))
        bigClean <- filter_poi(big, Classification != LASNOISE)
        writeLAS(bigClean, file = paste0(dirPath, groundPath, fileFinder, "_sor.laz"))
        safeBig <- big
        big <- bigClean
        rm(bigClean)
      }


      # DEM - 0.1 res and 0.01 disk size ####
      dem.res <- 0.1
      dem.p2r <- 0.01

      {
        cat(paste0("(*) CANOPY DEM res",dem.res*100,"cm discblur",dem.p2r*100,"cm... "))
        co <- capture.output(canop <- grid_canopy(big, res = dem.res, p2r(dem.p2r)))
        #plot(canop, main = paste0(groundCutHeight," DEM res",dem.res, " p2r",dem.p2r))

        if(rainFilter > 0){
          rf1 <- Sys.time()
          for(i in 1:rainFilter){
            if(i == 1){
              cat("RAIN FILTER ")
              normi <- big
            }
            cat(i)
            cat("x ")
            normi <- normalize_height(normi, canop)
            normi <- filter_poi(normi, Z <= -0.01)
            normi <- unnormalize_height(normi)
            co <- capture.output(canop <- grid_canopy(normi, res = dem.res, p2r(dem.p2r)))
            #plot(canop, main = paste0(groundCutHeight," DEM res",dem.res, " p2r",dem.p2r))
          }
          rf2 <- Sys.time()
          timeRainFilter <- as.difftime(rf2 - rf1)
          cat("DONE IN", round(timeRainFilter,1), units(timeRainFilter),"\n")
        }
      }
      writeRaster(canop, paste0(dirPath, groundPath, fileFinder, "_dem.tif"), overwrite = TRUE)
      writeRaster(canop, paste0(dirPath, groundPath, fileFinder, "_dem.grd"), overwrite = TRUE)

      # NDOM ####
      tryCatch(
        {

      cat(paste0("(*) NDOM in assembly... "))
      ndom <- canop - dtm_y
      ndom@extent
      #ndom@data@values[is.na(ndom@data@values)] <- 0
      ndom@data@values[!is.na(ndom@data@values) & ndom@data@values < 0] <- 0
      ndom@data@values[!is.na(ndom@data@values) & ndom@data@values > 44] <- NA
      ndom@data@values[!is.na(ndom@data@values) & ndom@data@values > 41] <- 41
      if(clip.radius != 0){
        #create a 40 m radius clipping
        cc_20 <- dismo::circles(data.frame(clip.x, clip.y), lonlat = F, clip.radius)
        ndom <- mask(crop(ndom, polygons(cc_20)),polygons(cc_20))
      }
      ndom@data@values[!is.na(ndom@data@values)][1:100] <- 41
      writeRaster(ndom, paste0(dirPath, groundPath, fileFinder, "_ndom.tif"), overwrite = TRUE)
      writeRaster(ndom, paste0(dirPath, groundPath, fileFinder, "_ndom.grd"), overwrite = TRUE)
      ndom@data@values[is.na(ndom@data@values)] <- 0


      png(filename = paste0(dirPath, imgPath, fileFinder,"_NDOM_raw.png"), 
          width = 800, height = 800, type = "cairo")
      par(xpd = F, mar = c(2,0,5,0), oma = c(0,0,0,0), xaxs='i', yaxs='i')
      if(clip.radius != 0){
        plot(0, type = "n",
             xlim = c(-clip.radius + clip.x, clip.radius + clip.x),
             ylim = c(-clip.radius + clip.y, clip.radius + clip.y), asp = 1, main = fileFinder, cex.main = 4, axes = F)
        raster::plot(ndom,  col = c(rainbow(19, start = 0.76, end = 0.60, rev = T), "#56008F", "#000000"), # every color is 1 m step (40 m span, 41 = taller than 40)
                     box = F, axes = F, legend = F,
                     asp = 1, add = T)
      } else {
        raster::plot(ndom,  col = c(rainbow(19, start = 0.76, end = 0.60, rev = T), "#56008F", "#000000"), # every color is 1 m step (40 m span, 41 = taller than 40)
                     box = F, axes = F, legend = F, cex.main = 4,
                     asp = 1, main = fileFinder)
      }
      axis(1, cex.axis = 2, lwd = 2)

      lgd_ <- c(">",seq(40, 1, by = -2))
      lgd_[nchar(lgd_) < 2] <- paste0("  ", lgd_[nchar(lgd_) < 2])
      legend("topright", title = "NDOM [m]",
             legend = lgd_, pt.cex = 3,
             fill = c("#000000", "#56008F", rainbow(19, start = 0.76, end = 0.60, rev = F)),
             border = NA,
             y.intersp = 0.7, cex = 1, text.font = 1, box.lty = 0)
      if(clip.radius != 0){
        clip(x1 = -clip.radius + clip.x + 0.5, x2 = clip.radius + clip.x - 0.5,
             y1 = -clip.radius + clip.y + 0.5, y2 = clip.radius + clip.y - 0.5)
      }
      abline(h = clip.y, v = clip.x, lwd = 1)
      if(clip.radius > 20){
        draw.circle(clip.x, clip.y, 20, lwd = 2, lty = 2)
      } else if(clip.radius > 10){
        draw.circle(clip.x, clip.y, 10, lwd = 2, lty = 2)
      }
      dev.off()


      }, error = function(error_condition) {
          warning("NDOM not created... Is that a problem?\n")
        })



      gstop <- Sys.time()
      cat("All DTMs done in a")
      print.difftime(round(gstop - gstart, 1))
    }



    ## writing laz files, sor filtering ####
    {
      cat("Writing out ground file ", paste0(dirPath, groundPath, fileFinder, "_ground.laz... \n"))
      writeLAS(ground, file = paste0(dirPath, groundPath, fileFinder, "_ground.laz"))
      LAS_ground <<- ground
      # groundSm <- decimate_points(ground, random(200))
      # groundSmPts <- groundSm@header@PHB$`Number of point records`
      # cat(" and", paste0(fileFinder, "_ground_sm.las... "))
      # writeLAS(groundSm, file = paste0(dirPath, groundPath, fileFinder, "_ground_sm.las"))
      #
      # rm(groundSm)
      if(!additionalSlices) rm(ground) # keep the ground for additional slices
      gc()


      vegetation <- filter_poi(big, Classification != 2) #0 = never classified, #1 = unclassified (vegetation) #18 = noise

      # "SOR" Noise filtering
      if(doFilter){
        cat("\nStarting with noise filtering (classify_noise with ivf algorithm in lidR):\n")
        cat("Voxel resolution:", filterRes, "\nminimum number of points in these voxels:", filterN, "\n")
        startSOR <- Sys.time()
        vegetation <- classify_noise(las = vegetation, algorithm = ivf(res = filterRes, n = filterN))
        stopSOR <- Sys.time()

        noiseP <- (sum(vegetation@data$Classification == 18))
        allP <- vegetation@header@PHB$`Number of point records`
        cat("Classified", thMk(noiseP), "points as noise, that is", round(noiseP/allP*100, 2), "%.\n")

        cat("Noise filtering done in a ")
        print.difftime(round(stopSOR - startSOR, 1))
        cat("\nWriting out")
      } else {
        cat("    and")
      }

      cat(" vegetation file", paste0(dirPath, groundPath, fileFinder, "_raw_veg.laz... "))
      writeLAS(vegetation, file = paste0(dirPath, groundPath, fileFinder, "_raw_veg.laz"))
      LAS_veg <<- vegetation
      LAS_veg_name <<- fileFinder
      cat("done!\n")


      {
        cat("\n The ground file up to",groundCutHeight*100,"cm contains", thMk(groundPts), "points, \n")
        cat("       the 5 cm voxelized file contains",thMk(groundFullVoxPts),"points, \n")
        cat("     and minimal voxel file up to 10 cm",thMk(groundMinVoxPts),"points.\n")
        cat("The raw vegetation file above contains", thMk(vegetation@header@PHB$`Number of point records`), "points.\n")
        cat("COMPARE: Comlete input cloud contains", thMk(big@header@PHB$`Number of point records`), "points (")
        cat(round(groundPts/big@header@PHB$`Number of point records`*100,0),"% ground)\n")
        if(1 == 2){
          #lidR::plot(vegetation)

          ints <- vegetation@data$Intensity
          cat("Intensity is scoping from",thMk(min(ints)),"to",thMk(max(ints)),"\n")
          #manual Histogram
          #hist(ints, freq = FALSE, breaks = c(0,(1:5)*10000,max(ints)))
          hist(ints, freq = FALSE)
          #plot(sort(ints[1:1000], decreasing = TRUE))
        }
      }

      rm(vegetation, ground, gr_vox)
      gc()

    }





    ## create handy slices ####
    if(exportSlice.upperLimit > 0){
      allStop <- Sys.time()
      cat("\nGround and vegetation splitting completed in ")
      print.difftime(round(allStop - allStart,1))

      # NORMALIZATION
      useCoarseGrid <- FALSE
      tryCatch(
        {
          cat("\nCreating additional slices - normalizing height... ")
          big <- normalize_height(big, dtm_y, na.rm = TRUE)
          cat("done!\n")
        }, error = function(error_condition) {
          cat("-> problem with the normalisation, the DTM model seems to be faulty...\n")
          useCoarseGrid <<- TRUE
        })

      if(useCoarseGrid){
        # take raster file 2, because first one produces some NAs at normalization
        cat("\"Roughly\" normalizing height (using the coarse 1 x 1 m dtm)...")

        tryCatch(
          {
            big <- normalize_height(big, dtm_rough, na.rm = TRUE) # need to save it in that intermediate object or it cannot unnormalize anymore
            cat("done!\n")
          }, error = function(error_condition) {
            cat("Error in creating the rough dtm model!\n")
            return()
          })
      }


      slicePath <- paste0(dirPath, groundPath, fileFinder,
                          "_clusterSlice_100to300.laz")
      #voxSlicePath <- paste0(dirPath, groundPath, fileFinder,
      #                       "_clusterSlice_100to300_vox.laz")
      slice <- filter_poi(big, Z < 3,
                          Z > 1)
      slice_un <- unnormalize_height(slice)

      cat("   o) slice for tree detection from 100 to 300 cm... ")
      writeLAS(slice_un, slicePath)
      cat("done!\n")

      if(exportSlice.lowerLimit != 1 | exportSlice.upperLimit != 3){

        slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",
                            exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,".laz")
        #voxSlicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",
        #                       exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,"_vox.laz")

        #rm(dtm, dtm_fine, dtm3)
        slice <- filter_poi(big, Z < exportSlice.upperLimit, Z > exportSlice.lowerLimit)
        slice_un <- unnormalize_height(slice)

        cat("   o) desired special slice ", exportSlice.lowerLimit*100, "-", exportSlice.upperLimit*100, " cm... ", sep = "")
        writeLAS(slice_un, slicePath)
        cat("done!\n")
      }





      if(additionalSlices){

        cat("   o) additional slices... 120-140... ")
        # 0 to 2
        slice <- filter_poi(big, Z < 1.4, Z > 1.2)
        slice_un <- unnormalize_height(slice)
        slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_120to140.laz")
        writeLAS(slice_un, slicePath)

        cat("50-200... ")
        # 0 to 2
        slice <- filter_poi(big, Z < 2, Z > 0.5)
        slice_un <- unnormalize_height(slice)
        slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_50to200.laz")
        writeLAS(slice_un, slicePath)

        cat("300-500... ")
        slice <- filter_poi(big, Z < 5, Z > 3)
        slice_un <- unnormalize_height(slice)
        slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_300to500.laz")
        writeLAS(slice_un, slicePath)


        if(groundMergeCut > 0){
          cutter <- (groundMergeCut/100)
          cat(paste0("0-", groundMergeCut, "... "))
          slice <- filter_poi(big, Z < cutter)
          slice_un <- unnormalize_height(slice)
          slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_0_ground_", groundMergeCut, ".laz")
          writeLAS(slice_un, slicePath)
        }

        cat("done!\n")

        slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,".laz")


      }
      rm(slice, slice_un)
    }
  }







  rm(big)
  gc()




  allStop <- Sys.time()
  cat("\nAll tasks completed. Global ")
  print.difftime(round(allStop - allStart,1))
  cat("\n\n")
  sink()

}











#' Takes the vegetation cloud from "extractVegetation" and seperates the stems
#' according to their intensity values
#'
#' @param fileFinder user defined name of the dataset
#' @param quantileIntensity in percent, how many most intense points are to be kept
#' @param CC_level a connected components parameter for fineness in cloudcompare (default: 10)
#' @param CC_numberPoints the minimum size in points of a component to not be dropped (default: 1000)
#' @param clipHeight in m if you need to set all stem seeds to a uniform height
#' @param bottomCut only in combination with clipHeight, in m lower border for cutting of relative heights, default to 0
#' @param bushPreparation if TRUE then there is additional density filtering to remove low branches and bushes
#' @param numberIntensity if you rather want to specify a number of most intense points to be kept
#' @param silent fewer outputs if set to TRUE
#' @param fast if TRUE, no detailled per stem information will be calculated
#'
#' @export
stemSplit <- function(fileFinder, quantileIntensity = 15, CC_level = 10, CC_numberPoints = 1000,
                      clipHeight = 3, bottomCut = 1, bushPreparation = FALSE, filterSOR = FALSE,
                      numberIntensity = 0, silent = TRUE, fast = TRUE, groundCutHeight = 0.5, amplitudeLowPass = 25.8,
                      cutWindow = c(-1000, -1000, 2000), retainPointClouds = FALSE, dirPath = paste0(getwd(),"/")){
  #### PREPARING SECTION ####

  start <- Sys.time()
  fileFinder <- removeUmlaut(fileFinder)

  dirPath <- paste0(getwd(),"/")
  metaInfo <- data.frame("infile" = fileFinder,
                         "qI" = 0,
                         "level" = CC_level,
                         "numP" = CC_numberPoints,
                         "minZ" = bottomCut,
                         "maxZ" = clipHeight,
                         "stems" = NA,
                         "maxPS" = NA,
                         "avPS" = NA,
                         "minPS" = NA,
                         "points" = NA,
                         "percSV" = NA,
                         "percSI" = NA,
                         "points.int" = NA,
                         "percIV" = NA,
                         "points.veg" = NA,
                         "qPoints" = 0,
                         "setString" = "",
                         "time" = 0)

  # Using CloudCompare to do Component Labelling
  level <-  CC_level
  numberOfPoints <-  CC_numberPoints

  ## EITHER FIXED NUMBER OF THRESHOLD INTENSITY, OR QUANTILE
  if(numberIntensity == 0){
    threshold_percent <- quantileIntensity
    metaInfo$qI <- quantileIntensity

    setString <- generateSetString(fileFinder, mode = "COMP", threshold_percent = threshold_percent,
                                   bottomCut = bottomCut, clipHeight = clipHeight,
                                   bushPreparation = bushPreparation, filterSOR = filterSOR,
                                   level = level, numberOfPoints = numberOfPoints, cutWindow = cutWindow,
                                   silent = silent)
    fileCode <- generateSetString(fileFinder, mode = "COMP", threshold_percent = threshold_percent,
                                  bottomCut = bottomCut, clipHeight = clipHeight,
                                  bushPreparation = bushPreparation, filterSOR = filterSOR,
                                  silent = silent)

  } else {
    threshold <- numberIntensity # getting all intensity values above 30.000
    metaInfo$qPoints <- numberIntensity

    setString <- generateSetString(fileFinder, mode = "COMP", threshold = threshold,
                                   bottomCut = bottomCut, clipHeight = clipHeight,
                                   bushPreparation = bushPreparation, filterSOR = filterSOR,
                                   level = level, numberOfPoints = numberOfPoints, cutWindow = cutWindow,
                                   silent = silent)
    fileCode <- generateSetString(fileFinder, mode = "COMP", threshold = threshold,
                                  bottomCut = bottomCut, clipHeight = clipHeight,
                                  bushPreparation = bushPreparation, filterSOR = filterSOR,
                                  silent = silent)
  }

  dbhPath <- paste0(dirPath,setString,"_dbh/")
  if(!dir.exists(dbhPath)) dir.create(dbhPath)

  sink(paste0(dbhPath,fileFinder,"_stemSplit_",format(Sys.time(), "%Y%m%d_%H%M"),"_Rcons.txt"), append = TRUE, split = TRUE)

  cat("Starting ELLFO Stem Component Detection for",fileFinder,"on",format(Sys.time()),"\n\n")
  #cat("\nLet's go with ALLGO cluster splitting :D!\n")
  cat("Working path is ",dbhPath,"\n")
  if(!exists("LAS_veg")){
    LAS_veg <<- NA
    LAS_veg_name <<- "blank"
  }


  {
    cat("Creating no folder structure yet, maybe later... done!\n")
    #
    # path.output.cluster.end <- paste0(dbhPath,"fineCluster/")
    # path.output.cluster.endgraph <- paste0(dbhPath,"graphSlice/")
    #
    # if(!dir.exists(paste0(dbhPath,"residuals/"))) dir.create(paste0(dbhPath,"residuals/"))
    # if(!dir.exists(path.output.cluster.end)) dir.create(path.output.cluster.end)
    # if(!dir.exists(path.output.cluster.endgraph)) dir.create(path.output.cluster.endgraph)
    #
    # if(referenced){
    #   plotReferenceFile(fileFinder = fileFinder, ref.file = ref, ref.plot_id = ref.plot_id,
    #                     cutWindow = cutWindow, writePNG = TRUE, pathPNG = dbhPath)
    # }
    # cat("done!\n")
  }




  options(warn=-1) #no warnings for you, Angela
  #But note that turning off warning messages globally might not be a good idea.
  #To turn warnings back on, use: options(warn=0)

  if(!silent){
    cat("     components are seperated in level",level,"\n")
    cat("     and we delete all components below",thMk(numberOfPoints),"points.\n")
    cat("     Floor starts at",groundCutHeight*100,"cm above ground.\n\n")
    cat("Or in a nutshell:",setString,"\n\n")

    cat("Can we start?\n")
    for(j in 1:5){
      Sys.sleep(1)
      cat("... ")
    }
    cat("\n\nLet's go!\n")
  } else {
    cat("Starting with:",setString,"\n\n")
  }
  metaInfo$setString <- setString



  if(2 == 1){



    if(file.exists(paste0(dbhPath,"slice_cluster.laz"))){
      cat("Skipping creation of cluster, loading old slice_cluster.laz file... ")
      co <- capture.output(sliVox <<- readLAS(paste0(dbhPath,"slice_cluster.laz")))
      cat("done!\n")
    } else {
      cat("Gathering diffuse slice file... \n")

      groundPath <- v.env$groundPath
      slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",bottomCut*100,"to",clipHeight*100,".laz")









    }
      # Creating the slice and the component file, see below

      #...

      ### VOXEL SYSTEMATICS #### NOT IMPLEMENTED YET
      vox.size <- 0.02
      cat("\nVoxelize points with a raster of",vox.size*100,"cm... ")
      {
        #t1 <- Sys.time()
        #thin4 <- voxelize_points(slice, vox.size)
        t2 <- Sys.time()
        #cat("lidR done")
        #print.difftime(t2-t1)
        thin5 <- TreeLS::tlsSample(slice, TreeLS::smp.voxelize(vox.size)) #bei TLS 0.015 # bei PLS ist 0.02 anders # auch da ist das mit smp.voxelize neu # iPAd 0.015
        t3 <- Sys.time()
        cat("done by treeLS. ")
        print.difftime(round(t3-t2, 1))
      }
      sliVox <<- thin5
      rm(slice, thin5)
      gc()






  }


  if(1 == 1){ ### old version
    groundPath <- v.env$groundPath
    slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",bottomCut*100,"to",clipHeight*100,".laz")
    slice <- NA

    #unlink(paste0(dirPath,setString), recursive = TRUE)
    #dir.create(dbhPath)
    #cat("new: ", fileFinder,"_intSeg_diffuse.las")
    #diffuseFile <- paste0(dirPath,substring(setString, 1, which(strsplit(setString, "_")[[2]]=="l")-2),"_intSeg_diffuse.las")
    diffuseFile <- paste0(dbhPath, fileCode,"_intSeg_diffuse.las")
    ClassifiedOutFileName <- paste0(setString,"_rawCloud_intensity.las") #only for CloudCompare
    dbhSidePath <- paste0(dbhPath, "cc_side/")
    dir.create(dbhSidePath)

    if(!file.exists(diffuseFile)){
      # NO DIFFUSE FILE, CHECK SLICE FILE FIRST

      ### CREATING SLICE FILE ####
      if(!file.exists(slicePath)){



        # CREATE RAW SLICE

        if(is.na(LAS_veg) || fileFinder != LAS_veg_name){
          cat("Reading in: ",paste0(dirPath,groundPath,fileFinder,"_raw_veg.las"),"...\n",sep = "")
          co <- capture.output(vegetation <- readLAS(file = paste0(dirPath,groundPath,fileFinder,"_raw_veg.las"), select = selector))

          if(retainPointClouds){
            cat("Retaining LAS_veg variable for ")
            LAS_veg <<- vegetation #retain big LAS file in memory, less performant
            LAS_veg_name <<- fileFinder
            cat(LAS_veg_name,"\n")
          }

        } else {
          cat("We use the prevailing \"LAS_veg\" variable to extract the",LAS_veg_name,"point cloud...\n")
          vegetation <- LAS_veg
        }



        tryCatch(
          {
            # read in raster file
            dtm1 <- raster(paste0(dirPath,groundPath,fileFinder,"_ground_min.grd"))
          }, error = function(error_condition) {
            cat("Error in reading the dtm-model, file not found!")
            return()
          })

        useCoarseGrid <- FALSE
        # NORMALIZATION
        tryCatch(
          {
            cat("Normalizing height... ")
            vegetation <- normalize_height(vegetation, dtm1, na.rm = TRUE)
            cat("done!\n")
          }, error = function(error_condition) {
            cat("-> problem with the normalisation, the DTM model seems to be faulty...\n")
            useCoarseGrid <- TRUE
          })

        if(useCoarseGrid){
          # read in raster file 2, because first one produces some NAs at normalization
          cat("\"Roughly\" normalizing height (using the coarse 1 x 1 m dtm)...")
          tryCatch(
            {
              # read in raster file
              dtm2 <- raster(paste0(dirPath,groundPath,fileFinder,"_ground_rough.grd"))
            }, error = function(error_condition) {
              cat("Error in reading the rough dtm-model, file not found!")
              return()
            })

          tryCatch(
            {
              vegetation <- normalize_height(vegetation, dtm2, na.rm = TRUE) # need to save it in that intermediate object or it cannot unnormalize anymore
              cat("done!\n")
            }, error = function(error_condition) {
              cat("Error in creating the rough dtm model!")
              return()
            })
        }

        #rm(dtm, dtm_fine, dtm3)
        slice <- filter_poi(vegetation, Z < clipHeight, Z > bottomCut)
        slice_un <- unnormalize_height(slice)
        cat("Creating global slice file at",slicePath,"... ")
        writeLAS(slice_un, slicePath)
        cat("done.\n")

        numP <- slice@header@PHB$`Number of point records`

        if(sum(cutWindow == c(-1000,-1000,2000)) == 3){
          cat("No cutting because of no cutWindow setting!\n")
        } else {
          cat("Just for output: clipping raw",thMk(numP),"points to dimensions", cutWindow, "leaves")
          XL <- cutWindow[1]
          YL <- cutWindow[2]
          width <- cutWindow[3]
          slice_un <- filter_poi(slice_un, X > XL, X < XL + width, Y > YL, Y < YL + width)
          numP <- slice_un@header@PHB$`Number of point records`
          cat(thMk(numP),"points.\n")
        }

        cat("Creating local cropped slice file at",slicePath,"... ")
        writeLAS(slice_un, paste0(dbhPath,"slice_raw.laz"))
        cat("done, retaining slice to do intensity component analysis.\n")

        slice <- slice_un
        rm(vegetation, slice_un)
        gc()

        stop <- Sys.time()
        print.difftime(round(stop - start, 1))


      } else # READING OLD SLICE FILE
      {
        cat("-> we use prevailing raw slice file, copying... ")
        if(sum(cutWindow == c(-1000,-1000,2000)) == 3){
          file.copy(slicePath, paste0(dbhPath,"slice_raw.laz"))
          cat("complete file copied.\n")
        } else {
          cat("reading in... ")
          co <- capture.output(slice <- readLAS(slicePath, select = selector))
          numP <- slice@header@PHB$`Number of point records`
          cat("found originally",thMk(numP), "points.\n")
          cat("Clipping to dimensions", cutWindow, "leaves ")

          XL <- cutWindow[1]
          YL <- cutWindow[2]
          width <- cutWindow[3]
          slice <- filter_poi(slice, X > XL, X < XL + width, Y > YL, Y < YL + width)

          numP <- slice@header@PHB$`Number of point records`
          cat(thMk(numP),"points, writing out slice_raw.laz... ")
          writeLAS(slice, paste0(dbhPath,"slice_raw.laz"))
          cat("done!\n")
        }
      }


      # CREATING DIFFUSE FILE
      if(is.na(slice)){
        cat("Reading in slice file... ")
        co <- capture.output(slice <- readLAS(slicePath, select = selector))
        cat("done, there were", thMk(slice@header@PHB$`Number of point records`), "points found.\n")
      }

      ### SOR FILTERING ###
      if(filterSOR){
        cat("Applying noise filter from inside point cloud (no separate settings specified):")
        slice <- filter_poi(slice, Classification < 2) #0 = never classified, #1 = unclassified (vegetation)  #18 = noise
        cat("remaining",thMk(slice@header@PHB$`Number of point records`),"points (approx.",round(slice@header@PHB$`Number of point records`/numP*100,1),"%).\n")
      }

      ### INTENSITY-FILTER ####
      cat("There were", thMk(slice@header@PHB$`Number of point records`), "points found.\n")
      metaInfo$points.veg <- slice@header@PHB$`Number of point records`
      cat("\nINTENSITY-FILTERING -")
      if(exists("threshold")){
        cat(" Absolute: \nIntensity higher than",thMk(threshold),"will be kept...\n")
      } else {
        threshold <- quantile(slice@data$Intensity,1 - threshold_percent/100) # all above are 5 %
        cat(" Relative: \nKeeping",round(threshold_percent,1),"% of all points with an intensity higher than",thMk(threshold),"...\n")
      }
      lasInt <- filter_poi(slice, Intensity > threshold)
      percentRem <- lasInt@header@PHB$`Number of point records`/slice@header@PHB$`Number of point records`
      metaInfo$points.int <- lasInt@header@PHB$`Number of point records`
      cat("There are", thMk(lasInt@header@PHB$`Number of point records`), "points remaining (equals only",
          round(percentRem*100,1),"% of original data).\n")

      if(lasInt@header@PHB$`X offset` < 0) lasInt@header@PHB$`X offset` <- 0 # correcting offsets, else cannot write las with negative offset... stupid rules...
      if(lasInt@header@PHB$`Y offset` < 0) lasInt@header@PHB$`Y offset` <- 0
      if(lasInt@header@PHB$`Z offset` < 0) lasInt@header@PHB$`Z offset` <- 0
      lasInt <- add_lasattribute(lasInt, 1:lasInt@header@PHB$`Number of point records`, "StemID", 'Stem ID Int Clust')
      #lidR::plot(lasInt)

      if(is.element("Amplitude", unlist(lasInt@header@VLR$Extra_Bytes$`Extra Bytes Description`))){
        lasInt <- filter_poi(lasInt, Amplitude < amplitudeLowPass)
        cat("Additional Amplitude filtering done!\n")
      }

      writeLAS(lasInt, file = diffuseFile)
      file.copy(diffuseFile, paste0(dbhSidePath,ClassifiedOutFileName))
      rm(threshold, slice)





    } else # READING IN OLD DIFFUSE FILE
    {
      cat("The diffuse file already exists (skipping intensity filtering part).\n")
      co <- capture.output(lasInt <- readLAS(file = diffuseFile, select = selector))
      cat("There were", thMk(lasInt@header@PHB$`Number of point records`), "points found.\n")
      file.copy(diffuseFile, paste0(dbhSidePath,ClassifiedOutFileName))
      cat("File",diffuseFile,"copied succesfully to the new directory.\n")
    }
  }











  #### PROCESSING CONNECTED COMPONENTS ####
  if(!silent)lidR::plot(lasInt, color = "Intensity")
  rm(lasInt)
  cat("\nUsing CloudCompare tool \"Label Connected Components\"...\n")
  if(fast){
    cat("Processing type is fast!\n")
  } else {
    cat("Processing type is slow but extensive!\n")
  }

  File.exe <- "Cloudcompare.exe"
  if(!exists("Path.CloudCompare")) Path.CloudCompare <- '"C:/Program Files/CloudCompare/"'
  Path.exe <- paste0(Path.CloudCompare,File.exe)



  #CloudCompare -O myhugecloud.bin -SS SPATIAL 0.1
  command.exe <- paste0(" ", Path.exe,
                        " -SILENT -o -GLOBAL_SHIFT AUTO ", dbhSidePath, ClassifiedOutFileName,
                        " -C_EXPORT_FMT LAS -EXTRACT_CC ",level," ",numberOfPoints
  )
  if(silent){
    system(command.exe, show.output.on.console = FALSE)
  } else {
    cat(command.exe,"\n")
    system(command.exe)
  }


  #### MERGING INTO SINGLE FILE #####

  # opening all files in the dir
  files <- list.files(path = dbhSidePath, full.names = TRUE)

  if(!file.exists(diffuseFile)){
    # we already did the filtering before, no need to create file again
    cat("Generating output file from",files[1],"\n")
    file.copy(files[1], diffuseFile)
  } else {
    #cat("No need to recreate output file\n")
  }
  file.remove(files[1]) # exclude input file in same folder
  files <- list.files(path = dbhSidePath, full.names = TRUE) # list again
  total <- length(files)
  metaInfo$stems <- total
  sizeList <- data.frame("stemID"=c(1:total))




  for(i in 1:total){
    cat("-",i,"of",total,"-")
    if(!silent)cat("read file",files[i])

    if(i %% 50 == 0){
      cat("Emptying mergedLAS...")
      if(!exists("mergedLASall")){
        mergedLASall <- mergedLAS
        rm(mergedLAS)
        cat(" new done!\n")
      } else {
        mergedLASall <- rbind(mergedLASall, mergedLAS)
        rm(mergedLAS)
        cat(" appended!\n")
      }
    }

    if(!exists("mergedLAS")){
      # First tree is exclusive, because we need to create new LAS object, change some day

      co <- capture.output(mergedLAS <- readLAS(file.path(files[i]), select = selector))
      mergedLAS$StemID <- i

      tempLAS <- mergedLAS
    } else {
      co <- capture.output(tempLAS <- readLAS(file.path(files[i]), select = selector))

      tempLAS@data$StemID <- i
      mergedLAS <- rbind(mergedLAS, tempLAS)
    }





    if(!fast){
      ## DISCOVER ALL POSSIBLE VARIABLES FROM SINGLE STEM FILE HERE
      heightSpan <- round(tempLAS@header@PHB$`Max Z` - tempLAS@header@PHB$`Min Z`,2)
      proj2d <- data.frame("x" = tempLAS@data$X, "y" = tempLAS@data$Y, "z" = tempLAS@data$Z)
      medx <- median(proj2d$x)
      medy <- median(proj2d$y)
      proj2d$dist <- apply(proj2d, 1, function(x) pointDistance(c(x[1],x[2]),c(medx, medy), lonlat = FALSE))
      #plot(sort(proj2d$dist, decreasing = TRUE)  )

      border <- quantile(proj2d$dist, 0.9) # removing 10 % of all furthest away points from center
      proj2d$remove <- proj2d$dist > border
      borderHeight <- round(tempLAS@header@PHB$`Min Z` + 2, 1) # only up to a height of 2 m
      medCenter <- proj2d[proj2d$remove == FALSE,] # here are only "valid" points, who are not far away
      medCenter <- medCenter[medCenter$z < borderHeight,] # and also lower than the border cutoff

      # coordinates of single stem
      sizeList$x[i] <- median(medCenter$x)
      sizeList$y[i] <- median(medCenter$y)
      sizeList$z[i] <- min(medCenter$z)

      #calculate distance for all points again from new absolute tree center
      medx <- median(medCenter$x)
      medy <- median(medCenter$y)
      proj2d$dist <- apply(proj2d, 1, function(x) pointDistance(c(x[1],x[2]),c(medx, medy), lonlat = FALSE))
      sizeList$minDist[i] <- quantile(proj2d$dist, 0.05) # getting the distance of 10 % closest points
      sizeList$maxDist[i] <- quantile(proj2d$dist, 0.95) # getting the distance of 10 % furthest points
      sizeList$heightSpan[i] <- heightSpan
      sizeList$numberOfPoints[i] <- tempLAS@header@PHB$`Number of point records`
    } else {
      #sizeList$x[i] <- NULL
      #sizeList$y[i] <- NULL
      #sizeList$z[i] <- NULL
      # coordinates of single stem
      sizeList$x[i] <- median(tempLAS@data$X)
      sizeList$y[i] <- median(tempLAS@data$Y)
      sizeList$z[i] <- median(tempLAS@data$Z)
      sizeList$minDist[i] <- NULL
      sizeList$maxDist[i] <- NULL
      sizeList$heightSpan[i] <- round(tempLAS@header@PHB$`Max Z` - tempLAS@header@PHB$`Min Z`,2)
      sizeList$numberOfPoints[i] <- tempLAS@header@PHB$`Number of point records`
    }

  }

  if(!exists("mergedLASall")){
    # less than 50 trees,
    mergedLASall <- mergedLAS
    rm(mergedLAS)
  } else {
    # more than 50, getting last trees into merged set
    if(exists("mergedLAS")){
      cat("Putting together last mergedLAS...")
      mergedLASall <- rbind(mergedLASall, mergedLAS)
      cat(" done!\n")
      rm(mergedLAS)
    }
  }


  #unique(mergedLAS@data$StemID)

  # Random Stem ID Part for viewing
  tryCatch(
    {
      mergedLASall <- add_lasattribute(mergedLASall, 0, "randomCol", "Random Stem ID")
      idlist <- unique(mergedLASall@data$StemID)
      rnlist <- sample(0:65535, length(idlist))
      mergedLASall@data$randomCol <- rnlist[match(mergedLASall@data$StemID, idlist)]
    }, error=function(e) cat("Something went wrong with random stem numbers...\n"))



  cat("\nCreated a stemset of",length(unique(mergedLASall@data$StemID)),"levels of trees.\n")
  tryCatch(writeLAS(mergedLASall, paste0(dbhPath,setString,"_intSeg_Stems.las")),
           error=function(e) cat("Something went wrong with outputting .las file...\n"))
  tryCatch(write.table(sizeList, paste0(dbhPath,setString,"_sizeList.txt"), row.names = FALSE, sep = "\t"),
           error=function(e) cat("Something went wrong with outputting .txt file...\n"))
  metaInfo$maxPS <- max(sizeList$numberOfPoints)
  metaInfo$minPS <- min(sizeList$numberOfPoints)
  metaInfo$avPS <- round(mean(sizeList$numberOfPoints),0)
  cat("Output file generation over.\n")
  lidR::plot(mergedLASall, color = "StemID")
  metaInfo$points <- mergedLASall@header@PHB$`Number of point records`
  metaInfo$percSV <- round(metaInfo$points / metaInfo$points.veg * 100,1)
  metaInfo$percSI <- round(metaInfo$points / metaInfo$points.int * 100,1)
  metaInfo$percIV <- round(metaInfo$points.int / metaInfo$points.veg * 100,1)
  metaInfo$time <- format(Sys.time(), "%y%m%d_%X")
  rm(mergedLASall)


  #### OUTPUT META-INFO ######


  cat("Writing output metaInfo...\n")
  tryCatch(
    {
      metaOld <- read.table(paste0(dirPath,"metaStemTrys.txt"), header = TRUE)
      #meta <- merge.data.frame(x = meta, y = metaOld, by = c("file","qI","qPoints"), all.x = TRUE)

      if(is.na(metaInfo$points.int)){
        location <- match(paste0(metaInfo$infile, metaInfo$qI, metaInfo$qPoints),
                          paste0(metaOld$infile, metaOld$qI, metaOld$qPoints))
        if(!is.na(location)) metaInfo$points.int <- metaOld$points.int[location]
      }
      if(is.na(metaInfo$points.veg)){
        cat("Missing veg points...\n")
        location <- match(paste0(metaInfo$infile, metaInfo$qI, metaInfo$qPoints),
                          paste0(metaOld$infile, metaOld$qI, metaOld$qPoints))
        if(!is.na(location)) metaInfo$points.veg <- metaOld$points.veg[location]
      }

      # S... stems (seperated with missing noises < numP)
      # I... intensity (thinned vegetation to qI)
      # V... vegetation (whole point cloud in)

      #for example percSV means what share of total vegetation (all points) is in the stem file

      metaInfo$percSV <- round(metaInfo$points / metaInfo$points.veg * 100,1)
      metaInfo$percSI <- round(metaInfo$points / metaInfo$points.int * 100,1)
      metaInfo$percIV <- round(metaInfo$points.int / metaInfo$points.veg * 100,1)
      indexDoubled <- match(metaInfo$setString, metaOld$setString)
      if(!is.na(indexDoubled)){
        cat("Already did this set, rewrite old one...\n")
        metaOld <- metaOld[-indexDoubled,]
      }
      meta <- rbind.data.frame(metaOld, metaInfo)
      meta <- meta[order(meta$numP),]
      meta <- meta[order(meta$level),]
      meta <- meta[order(meta$qI),]
      meta <- meta[order(meta$infile),]
      write.table(meta, paste0(dirPath,"metaStemTrys.txt"), row.names = FALSE, sep = "\t")


    },
    error=function(e){
      if(file.exists(paste0(dirPath,"metaStemTrys.txt"))){
        write.table(metaInfo,paste0(dirPath,"metaStemTrys_ERR.txt"), row.names = FALSE, sep = "\t")
        cat("Error in metadata writing.")
      } else {
        cat("Need to create new metadata-file.\n")
        write.table(metaInfo,paste0(dirPath,"metaStemTrys.txt"), row.names = FALSE, sep = "\t")
      }
    }
  )

  cat("Done.\n")



  if(!retainPointClouds){
    LAS_veg <<- NA #removing old LAS_veg to ensure performance
    gc()
  }

  cat("\nEnd of routine.\n")
  stop <- Sys.time()
  print.difftime(round(stop - start,1))
  cat("Deleting conComp folder...", dbhSidePath,"\n")
  unlink(substr(dbhSidePath,1,nchar(dbhSidePath)-1), recursive = TRUE)
  rm(start, stop)
  #print(metaInfo)
  cat("Case",setString,"is closed.\n\n")
  options(warn=0)

  sink()
}



