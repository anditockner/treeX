
if(!exists("LAS_veg")){
  LAS_veg <<- NA
  LAS_ground <<- NA
  LAS_veg_name <<- "blank"
  }





#' Ground classification of raw point cloud
#'
#' Reads in a .las or .laz file from disk and extracts the vegetation cloud
#' using the cloth simulation filter (scf) via lasground from the lidR package
#' NEW: Can create an initial slice already! for the first faster clustering step
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
                              doFilter = FALSE, filterRes = 0.05, filterN = 27, tooBig = F,
                              clip.radius = 0, clip.trajectory.distance = 0,
                              clip.x = 0, clip.y = 0,
                              selector = "xyzcit0", dtm.path = "",
                              trafoMatrix.path = "",
                              rainFilter = 0,
                              additionalSlices = TRUE,
                              exportSlice.upperLimit = 3, exportSlice.lowerLimit = 1, #m
                              exportClippedLAS = FALSE,
                              dirPath = paste0(getwd(), "/")){


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

    tooBig = F
    groundPath <- paste0("_total_ground_veg/")
  }


  groundPath <- v.env$groundPath
  if(!dir.exists(paste0(dirPath, groundPath))) dir.create(paste0(dirPath, groundPath))
  imgPath <- "png/"
  if(!dir.exists(paste0(dirPath, groundPath, imgPath))) dir.create(paste0(dirPath, groundPath, imgPath))


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


    # slicePath <- paste0(dirPath, "slice_", exportSlice.lowerLimit,
    #           "to", exportSlice.upperLimit, "/")
    slicePath <- groundPath



    if(!dir.exists(slicePath)) dir.create(slicePath)
    # if(!dir.exists(paste0(out.path20, groundPath))) dir.create(paste0(out.path20, groundPath))
    # if(!dir.exists(paste0(out.path40, groundPath))) dir.create(paste0(out.path40, groundPath))
    #
    # inputFile.list <- list.files(inputFile.path, pattern = ".laz")
    # inputFile <- "none.laz"

  }


  sink(paste0(dirPath, groundPath, fileFinder, "_extractVegetation_", format(Sys.time(), "%Y%m%d_%H%M"), "_Rcons.txt"), append = TRUE, split = TRUE)
  cat("Starting to create a ground and vegetation model in set", fileFinder, "\n")
  cat("Today is", format(Sys.time()), "\n")
  if(clip.trajectory.distance){
    cat(paste0("Clipping input file to trajectory +", clip.trajectory.distance, "m.\n"))
  }
  cat("Working at", dirPath, "\n\n")
  #if(!dir.exists(dirPath)) dir.create(dirPath)




  allStart <- Sys.time()

  if(tooBig){
    cat("tooBig: Reading only every 3rd point: ", LASfile, "...\n", sep = "")
    big <- readLAS(paste0(LASfile), select = selector, filter = "-keep_every_nth 3")
    if(clothSize == 0.1){
      clothSize <- 0.5
    }
  } else {
    cat("Reading in: ", basename(LASfile), "... \n", sep = "")
    cat("(from ", dirname(LASfile), ")...", sep = "")
    if(clip.radius > 0){
      cat("\nGoing to produce a circle of r =",clip.radius,"m")
      if(clip.x != 0 || clip.y != 0){
        cat(paste0(" - base point is set to x=", round(clip.x, 2), " and y=", round(clip.y, 2), "m"))
      }
      cat("... ")
      big <- readLAS(paste0(LASfile), select = selector, filter = paste0("-keep_circle ", clip.x, " ", clip.y, " ",clip.radius))
    } else {
      big <- readLAS(paste0(LASfile), select = selector)
    }
    cat("done!\n")
  }


  if(!is.element("Intensity", colnames(big@data))){
    big@data$Intensity <- 5L
    warning("NO INTENSITY VALUES PRESENT!!!\n\n")
  }
  cat("Checking if intensity values are present...\n")
  if(max(big@data$Intensity) == 0){
    big@data$Intensity <- 5L
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



  # copy trajectory
  oneDone <- FALSE
  txtExists <- FALSE
  plyExists <- FALSE
  tryCatch({
    tempName <- basename(LASfile)

    finLaz <- strfind(tempName, "_100pct")
    if(is.null(finLaz)){
      finLaz <- strfind(tempName, ".las")
    }
    if(is.null(finLaz)){
      finLaz <- strfind(tempName, ".laz")
    }

    
    
    
    cat("TRAJECTORY: Trying to find file in",paste0(dirname(LASfile), "/", substr(tempName, 1, finLaz-1),"...\n    ..._results_traj.txt "))
    
    
    trajfile <- paste0(dirname(LASfile), 
                       "/", substr(tempName, 1, finLaz-1), "_results_traj.txt")
    if(!file.exists(trajfile)){
      trajfile <- paste0(dirname(LASfile), 
                         "/", substr(tempName, 1, finLaz-1), ".gs-traj")
    }
    if(!file.exists(trajfile)){
      trajfile <- paste0(dirname(LASfile), 
                         "/", substr(tempName, 1, finLaz-1), ".txt")
    }
    if(!file.exists(trajfile)){
      trajfile <- paste0(dirname(LASfile), 
                         "/", substr(tempName, 1, finLaz-1), "_traj.txt")
    }
    if(!file.exists(trajfile) & clip.trajectory.distance > 0){
      return("No trajectory found! Clipping cannot be done...\n\n")
    }
    
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

  }, error = function(error_condition) {
    cat("-> problem with the trajectory! ")
    if(oneDone) cat("Only the .txt file was copied!\n")
  })

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
      cat("")
    }, error = function(error_condition) {
      cat("NO TRANSFORMATION WAS DONE!!!\n\n\n")
      warning("Transformation process was not sucessful!")
    })
  }


  cat("\n")




  # DTM SECTION ####

  #dtm.path <- "D:/_later/Stefl_need_automatic_coregistration/_total_ground_veg/JA2_ground_min.grd"
  #big <- readLAS("D:/_later/Stefl_need_automatic_coregistration/_total_ground_veg/JA2_clusterSlice_50to200.laz", select = "xyzcit0")
  #big <- decimate_points(big, random(10))


  
  if(clip.trajectory.distance > 0 && txtExists){
    clipTime <- Sys.time()
    cat("Reading trajectory for clipping input... ")
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
    
    


    traj2 <- traj[seq(from = 1, to = nrow(traj), by = 50),]


    # trajectory creates usually 83 points per second,
    # we reduce them to 2 points per second for alpha-hulling
    #traj$X.time[2] - traj$X.time[1] original spacing of points is 1/100 of a second
    #traj2$X.time[2] - traj2$X.time[1] # now trajectory is spaced every 0.5 seconds
    traj <- traj2
    traj$col <- rainbow(length(traj[,1]), end = 0.7, rev = T)
    rm(traj2)

    #dups <- duplicated.data.frame(data.frame("x" = traj$x, "y" = traj$y))
    shiftX <- min(traj$x)
    shiftY <- min(traj$y)
    if(shiftX > 1000){
      cat("Shifting x by", shiftX, "m...\n")
      traj$x <- traj$x - shiftX
    }
    if(shiftY > 1000){
      cat("Shifting y by", shiftY, "m...\n")
      traj$y <- traj$y - shiftY
    }


    cat("Clipping total cloud to trajectory +", clip.trajectory.distance, "m radius (a=25m)...")
    borderHull <- ashape(traj$x, traj$y, alpha = 25)
    owin.window <-owin(xrange=range(traj$x), yrange=range(traj$y))
    trees_edges_ppp <- psp(borderHull$edges[,3], borderHull$edges[,4], borderHull$edges[,5], borderHull$edges[,6], window=owin.window)
    hull.traj <- dilation(trees_edges_ppp, r=clip.trajectory.distance)
    owin.traj <- owin(poly=hull.traj$bdry[[1]])
    area.traj <- area.owin(owin.traj)

    big_sm <- decimate_points(big, random(1))
    #big_sm <- readLAS("D:/Ebensee/_in_laz/097_20065_1911_0730_100pct_scan.laz", select = "xyzcit0")
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

#
# plot(outer)
#     #big <- big_sm
# big@data$X <- big@data$X + shiftX
# big@data$Y <- big@data$Y + shiftY
    pointsBefore <- big@header@PHB$`Number of point records`
    big <- clip_roi(big, outer)
    pointsAfter <- big@header@PHB$`Number of point records`

    pointsLost <- pointsBefore - pointsAfter
    cat("done!\n")
    cat("Remain", thMk(pointsAfter), "pts after clipping to shape.\n")
    cat("We lost", thMk(pointsLost), "pts (or", round(pointsLost/pointsBefore*100,1), "% of original pts)\n")

    png(paste0(dirPath, groundPath, imgPath, fileFinder, "_traj_clipping.png"),
        height = diff(hull.traj$yrange)*5, width = diff(hull.traj$xrange)*5)
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


    if(exportClippedLAS){
      writeLAS(big, paste0(out.path.clip, fileFinder, "_traj", clip.trajectory.distance, "m.laz"))
    }
    gc()
  }


  # DTM GIVEN ####
  if(file.exists(dtm.path)){
    # if DTM is given, we only need to cut it once to absolute values (need no puffer)

    cat("\nDTM-Section:\n")
    gstart <- Sys.time()
    tryCatch(
      {
        # read in raster file
        cat("Reading input dtm from",dtm.path,"\n")
        dtm_y <- raster(dtm.path)
      }, error = function(error_condition) {
        warning("Error in reading the dtm-model, file", basename(dtm.path), "not found!")
        return()
      })
    dtCrop <- crop(dtm_y, extent(c(big@header$`Min X`, big@header$`Max X`, big@header$`Min Y`, big@header$`Max Y`)))



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

    png(filename = paste0(dirPath, groundPath, imgPath, fileFinder,"_NDOM_raw.png"), width = 800, height = 800)
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
    if(clip.radius > 10){
      draw.circle(clip.x, clip.y, 10, lwd = 2, lty = 2)
    } else if(clip.radius > 20){
      draw.circle(clip.x, clip.y, 20, lwd = 2, lty = 2)
    }
    dev.off()


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
      gr_vox <- tlsSample(ground, smp.voxelize(vox.gr))
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


    # create handy slices ####
    if(exportSlice.upperLimit > 0){
      allStop <- Sys.time()
      cat("\nGround and vegetation splitting completed in ")
      print.difftime(round(allStop - allStart,1))
      cat("\nCreating a slice from", exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,"cm: \n")

      # NORMALIZATION IS ALREADY DONE!
      slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",
                          exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,".laz")
      voxSlicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",
                             exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,"_vox.laz")

      #rm(dtm, dtm_fine, dtm3)
      slice <- filter_poi(big, Z < exportSlice.upperLimit, Z > exportSlice.lowerLimit)
      slice_un <- unnormalize_height(slice)
      cat("Creating global slice file at", slicePath, "... ")
      writeLAS(slice_un, slicePath)
      cat("done.\n")

      if(additionalSlices){
        cat("Creating additional slices... 120-140... ")
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
    ## old stuff - creating new ground model and so on...
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
          gr_vox <- tlsSample(ground, smp.voxelize(vox.gr))
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
        gr_vox <- tlsSample(ground, smp.voxelize(vox.gr))
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
        gr_vox <- tlsSample(ground, smp.voxelize(vox.gr))
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
        # would be nice, but is faaaar to slow (takes up to 10 hours for 0.5 ha)
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
      #big <- readLAS("F:/_Jauch_ALL/in_sf/Jauch/0_Bach_Boendl/2023-08-16_15-22-11_Jauch_bei_Bach_regen_100pct_scan.laz", select = "xyzit")

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
      #writeRaster(canop, paste0("D:/canop_final.tif"), overwrite = TRUE)
      }
      writeRaster(canop, paste0(dirPath, groundPath, fileFinder, "_dem.tif"), overwrite = TRUE)
      writeRaster(canop, paste0(dirPath, groundPath, fileFinder, "_dem.grd"), overwrite = TRUE)

      # NDOM ####
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


      png(filename = paste0(groundPath, imgPath, fileFinder,"_NDOM_raw.png"), width = 800, height = 800)
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
      if(clip.radius > 10){
        draw.circle(clip.x, clip.y, 10, lwd = 2, lty = 2)
      } else if(clip.radius > 20){
        draw.circle(clip.x, clip.y, 20, lwd = 2, lty = 2)
      }
      dev.off()




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


      vegetation <- filter_poi(big, Classification < 2) #0 = never classified, #1 = unclassified (vegetation) #18 = noise

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
      cat("\nCreating a slice from", exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,"cm: \n")

      # NORMALIZATION
      useCoarseGrid <- FALSE
      tryCatch(
        {
          cat("Normalizing height... ")
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

      # previous: load in old DTMs, deprecated
      {
        # if(!noDTMs){
        #   # all that above
        # } else {
        #   # load external DTMs
        #   cat("Loading previous dtms...\n")
        #   tryCatch(
        #     {
        #       cat("DTM Using input model from", paste0(dirPath,groundPath,fileFinder,"_ground_c1.grd...\n"))
        #       # read in raster file
        #       dtm_fine <- raster(paste0(dirPath,groundPath,fileFinder,"_ground_c1.grd"))
        #     }, error = function(error_condition) {
        #       cat("Error in reading the dtm-model, try the second one...")
        #       useCoarseGrid <<- TRUE
        #     })
        #   if(!useCoarseGrid){
        #     tryCatch(
        #       {
        #         cat("Normalizing height... ")
        #         vegetation <- normalize_height(vegetation, dtm_fine, na.rm = TRUE)
        #         cat("done!\n")
        #       }, error = function(error_condition) {
        #         cat("-> problem with the normalisation, the DTM model seems to be faulty...\n")
        #         useCoarseGrid <<- TRUE
        #       })
        #   }
        #
        #   if(useCoarseGrid){
        #     # take raster file 2, because first is not available or produces some NAs at normalization
        #     cat("DTM Using rough 1 x 1 m input model from", paste0(dirPath,groundPath,fileFinder,"_ground_a2.grd...\n"))
        #     tryCatch(
        #       {
        #         # read in raster file
        #         dtm_rough <- raster(paste0(dirPath,groundPath,fileFinder,"_ground_a2.grd"))
        #       }, error = function(error_condition) {
        #         cat("Error in reading the dtm-model, file not found!")
        #         cat("\n\n--> There will be warnings and the slices will be mishapped!!!\n\n")
        #       })
        #
        #     tryCatch(
        #       {
        #         cat("Normalizing roughly height... ")
        #         vegetation <- normalize_height(vegetation, dtm_rough, na.rm = TRUE) # need to save it in that intermediate object or it cannot unnormalize anymore
        #         cat("done!\n")
        #       }, error = function(error_condition) {
        #         cat("Error in creating the rough dtm model!\n")
        #       })
        #   }
        # }
      }


      slicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",
                          exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,".laz")
      voxSlicePath <- paste0(dirPath,groundPath,fileFinder,"_clusterSlice_",
                             exportSlice.lowerLimit*100,"to",exportSlice.upperLimit*100,"_vox.laz")


      #rm(dtm, dtm_fine, dtm3)
      slice <- filter_poi(big, Z < exportSlice.upperLimit, Z > exportSlice.lowerLimit)
      slice_un <- unnormalize_height(slice)
      cat("Creating global slice file at", slicePath, "... ")
      writeLAS(slice_un, slicePath)
      cat("done.\n")

      if(additionalSlices){

        cat("Creating additional slices... 120-140... ")
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

