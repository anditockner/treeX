globalTimeDiff <<- "no timez"



#' @export
changeLASVEG <- function() {
  groundPath <- v.env$groundPath
  cat(groundPath, " is ground path\n")
  cat(LAS_veg_name, " is LAS veg name")
  cat(outLAS_name, " is LAS out name")
  # LAS_veg <<- "hi"
}


#' Tree crown segmentation
#'
#' Classifies all trees and isolates the crown,
#' stem coordinates are required from trees_dbh.txt
#' which are used as seed points, from there it grows the stems
#' until all points are allocated to a stemID
#'
#'
#' @param fileFinder user defined name of the dataset
#' @param cutWindow c(xL, yL, width) square frame to detect trees from
#' @param doReferencedOnly deletes all the above 20000 stems, they do not make senseful trees
#' @param zScale relative number of Z-variable multiplication for a better height progress
#' @param merged set TRUE if using a combined stem-file with field "size" for combining
#' @param totalRuns number of iterations for discovering crowns, default 500
#' @param limitShare if we are adding less points from blanks than that, then incrementing searching distance (default: 0.005 = 0.5 %)
#' @param zScale relative number of Z-variable multiplication for a better height progress
#' @param voxelSize before the region growing the point cloud can be voxelized, in cm, default = 0 (no voxelisation)
#' @param incrementDistance the step of raising searching distance if less than limitShare points are found
#'
#' @param CC_level a connected components parameter for fineness in cloudcompare (default: 10)
#' @param CC_numberPoints the minimum size in points of a component to not be dropped (default: 1000)
#' @param clipHeight in m if you need to set all stem seeds to a uniform height
#' @param bottomCut only in combination with clipHeight, in m lower border for cutting of relative heights, default to 0
#' @param bushPreparation if TRUE then there is additional density filtering to remove low branches and bushes
#' @param filterSOR if TRUE then we use the SOR clean file
#' @param numberIntensity if you rather want to specify a number of most intense points to be kept
#' @param silent fewer outputs or plots if set to TRUE
#' @param useTreeFile loading and external trees_dbh.txt file
#' @param fast if TRUE, no detailled per stem information will be calculated
#' @param retainPointClouds if TRUE, the used point clouds will be stored for another run (more performant for different settings), if FALSE the point clouds will be passed to next function layer (always) and then discarded.
#'
#' @export
crownFeel <- function(fileFinder, cutWindow = c(-1000, -1000, 2000), ipad = FALSE, doReferencedOnly = FALSE, referenced = FALSE,
                      zScale = 2, merged = FALSE, totalRuns = 1000, limitShare = 0.004, incrementDistance = 0.005,
                      limitStems = 0, tileClipping = 5, diagonals = FALSE,
                      useTreeFile = "",
                      selector = "xyzcit0", distanceCounter_cm_limit = 220,
                      quantileIntensity = 15, CC_level = 10, CC_numberPoints = 1000, durationMins = 2400, maximumDistance = 0.50,
                      clipHeight = 3, bottomCut = 1, bushPreparation = FALSE, filterSOR = FALSE, voxelSize = 0,
                      numberIntensity = 0, silent = TRUE, fast = TRUE, plotEvery = 100, mode = "ALLGO",
                      frame.up = 0.3, frame.down = 0.3, frame.rad = 1.05, # cutting m up and down for start seed, frame.rad factor 1.6xDBH
                      retainPointClouds = FALSE, dirPath = paste0(getwd(), "/")) {
  if (ipad) {
    clipHeight <- 2.2
    bottomCut <- 0.2
    selector <- "xyzcitRBG0"
  }

  # Creates values with thousandmark (are prettier to look at)
  # thMk <- function(val) {
  # val2 <- format(val, big.mark = ".", decimal.mark = ",", scientific = FALSE)
  # return(val2)
  # }


  ## TROUBLESHOOTING #
  if (2 == 1) {
    diagonals <- FALSE
    distanceCounter_cm_limit <- 220
    selector <- "xyzcit0"
    useTreeFile <- ""
    cutWindow <- c(-1000, -1000, 2000)
    tileClipping <- TRUE
    zScale <- 2
    merged <- FALSE
    totalRuns <- 1000
    limitShare <- 0.004
    limitStems <- 0
    incrementDistance <- 0.005


    quantileIntensity <- 15
    CC_level <- 10
    CC_numberPoints <- 1000
    durationMins <- 1200
    maximumDistance <- 0.50

    clipHeight <- 3
    bottomCut <- 1
    bushPreparation <- FALSE
    filterSOR <- FALSE
    voxelSize <- 10

    numberIntensity <- 0
    silent <- TRUE
    fast <- TRUE
    plotEvery <- 100
    mode <- "ALLGO"
    retainPointClouds <- FALSE
  }

  # PREPARATION SECTION
  {
    fileFinder <- removeUmlaut(fileFinder)

    if (!exists("LAS_veg")) {
      LAS_veg <<- NA
      LAS_veg_name <<- "blank"
    }

    if (!exists("outLAS")) {
      outLAS <<- NA # all crowns together in one las file after crownFeeling()
      outLAS_name <<- ""
    }

    if (!exists("blankLAS")) blankLAS <<- NA

    groundPath <- "_total_ground_veg/"
    # groundPath <- v.env$groundPath
    gstart <- Sys.time()
    mergedSet <- FALSE



    howManyClosestPoints <- 1 # finding 1 nearest neighbour to every blank point
    # CAUTION: changing it to more than 1 is highly inefficient...


    locationStr <- generateSetString(cutWindow = cutWindow, silent = silent)
    XL <- cutWindow[1]
    YL <- cutWindow[2]
    width <- cutWindow[3]


    if (mode == "COMP") {
      setString <- generateSetString(fileFinder,
        mode = mode,
        threshold = numberIntensity, threshold_percent = quantileIntensity,
        bottomCut = bottomCut, clipHeight = clipHeight,
        bushPreparation = bushPreparation, filterSOR = filterSOR,
        level = CC_level, numberOfPoints = CC_numberPoints, silent = silent
      )
      dbhPath <- paste0(dirPath, setString, "_dbh/")

      if (!dir.exists(dbhPath)) {
        setString <- generateSetString(fileFinder,
          mode = mode, cutWindow = cutWindow,
          threshold = numberIntensity, threshold_percent = quantileIntensity,
          bottomCut = bottomCut, clipHeight = clipHeight,
          bushPreparation = bushPreparation, filterSOR = filterSOR,
          level = CC_level, numberOfPoints = CC_numberPoints, silent = silent
        )
        dbhPath <- paste0(dirPath, setString, "_dbh/")
      }
    } else if (mode == "ALLGO") {
      dbhPath <- generateSetString(
        fileFinder = fileFinder, mode = mode,
        clipHeight = clipHeight, bottomCut = bottomCut,
        bushPreparation = bushPreparation,
        filterSOR = filterSOR, silent = TRUE
      )
      dbhPath <- paste0(dirPath, dbhPath, "_dbh/")

      if (!dir.exists(dbhPath)) {
        dbhPath <- generateSetString(
          fileFinder = fileFinder, mode = mode, cutWindow = cutWindow,
          clipHeight = clipHeight, bottomCut = bottomCut,
          bushPreparation = bushPreparation,
          filterSOR = filterSOR, silent = TRUE
        )
        dbhPath <- paste0(dirPath, dbhPath, "_dbh/")
      }

      if (useTreeFile != "") {
        if (!dir.exists(dbhPath)) {
          cat("Creating directory", dbhPath, "\n")
          dir.create(dbhPath, recursive = T)
        }
        if (!file.exists(useTreeFile)) {
          warning("There is an input file specified that could not be found!")
        } else {
          file.copy(
            from = useTreeFile,
            to = paste0(dbhPath, "trees_dbh.txt"), overwrite = T
          )
        }
      }

      if (!dir.exists(dbhPath)) {
        warning("There are no clustered input files! \nPlease run functions extractVegetation and clusterSplit / stemSplit before!\n")
        cat("Terminating on", format(Sys.time()), "\n\n")
      }
    }
    crownString <- generateSetString(limitShare = limitShare, zScale = zScale, voxelSize = voxelSize)
    crownPath <- paste0(dbhPath, "crowns", crownString, "/")

    if (!dir.exists(crownPath)) dir.create(crownPath)

    globalStart <- Sys.time()
    startTimeString <- format(globalStart, "%Y%m%d_%H%M")
    sink(paste0(crownPath, fileFinder, "_crownFeel_", startTimeString, "_Rcons.txt"), append = TRUE, split = TRUE)

    cat("Starting Crown Segmentation for", fileFinder, "on", format(Sys.time()), "\n\n")


    if (exists("dbhPath")) cat("dbhPath generated: ", dbhPath, "\n")
    cat("crownString generated: ", crownString, "\n")
    cat("locationString generated:", locationStr, "\n\n")
    # cat("limitShare:", limitShare, "\n") # how much percent of points need to be added per round
    # if actual points added lower than limitShare, then searchDistance is increased to find more points
    # limitShare too low (0.0001 %) will lead to a low searchDistance (mms)
    # and for too sparse point clouds there will be no good allocations, as
    # a dense tree will absorb all other trees
    # limitShare too high (0.01 %) will lead to a high searchDistance and to a merging between neighbors
  }

  if (is.na(LAS_veg) || LAS_veg_name != fileFinder) {
    totalCloud.name <- paste0(dirPath, groundPath, fileFinder, "_raw_veg.laz")
    groundCloud.name <- paste0(dirPath, groundPath, fileFinder, "_ground.laz")
    cat("Reading in totalCloud from", totalCloud.name, "\n")
    totalCloud <- readLAS(totalCloud.name, select = selector)
    cat("Reading in groundCloud from", groundCloud.name, "\n")
    groundCloud <- readLAS(groundCloud.name, select = selector)

    if (retainPointClouds) {
      cat("Retaining LAS_veg and LAS_ground variable for ")
      LAS_veg <<- totalCloud # retain big LAS file in memory, less performant
      LAS_ground <<- groundCloud
      LAS_veg_name <<- fileFinder
      cat(LAS_veg_name, "\n")
    }
  } else {
    cat("We use the prevailing \"LAS_veg\" variable to extract the", LAS_veg_name, "point clouds...\n")
    totalCloud <- LAS_veg
    if (is.na(LAS_ground)) {
      cat("Reading in groundCloud from", groundCloud.name, "\n")
      groundCloud <- readLAS(groundCloud.name, select = selector)
    } else {
      cat("Also using old groundCloud!\n")
      groundCloud <- LAS_ground
    }
  }

  vegCloud <- totalCloud
  totalCloud <- rbind(totalCloud, groundCloud)


  if (retainPointClouds == FALSE) {
    LAS_veg <<- NA # removing old LAS_veg to ensure performance
    gc()
  }


  if (mode == "COMP") {
    #### READING IN DATA ####
    stemCloud.name <- paste0(dbhPath, setString, "_intSeg_Trees.las")
    if (!file.exists(stemCloud.name)) {
      stemCloud.name <- paste0(dbhPath, setString, "_intSeg_Stems.las")
      cat("\nThere was no _intSeg_Trees.las found! No referencing added.\n")
      cat("Crown feeling is done only with stem information.\n")
    }
    cat("Reading in stemCloud from", stemCloud.name, "\n")
    stemCloud <- readLAS(stemCloud.name, select = selector)


    idlist <- unique(stemCloud@data$StemID)
    try(rnlist <- stemCloud@data$randomCol[match(idlist, stemCloud@data$StemID)])

    totalCloud <- add_lasattribute_manual(totalCloud, 0, "StemID", "Single Stem ID", type = "short")
    if (exists("rnlist")) totalCloud <- add_lasattribute_manual(totalCloud, 0, "randomCol", "Random Stem ID", type = "short")
    if (mergedSet) totalCloud <- add_lasattribute(totalCloud, 0, "size", "Cluster size")
    # lidR::plot(stemCloud)

    # length(stemCloud@data$X)
    cat("Merging...\n")
    totalCloud@data <- rbind(totalCloud@data, stemCloud@data, use.names = TRUE)
    stemCloud <- totalCloud
    stemCloud <- add_lasattribute_manual(stemCloud, 999, "runJay", "number of run added to total", type = "short")

    rm(totalCloud)
    gc()
  } else if (mode == "ALLGO") {
    clustList.name <- paste0(dbhPath, "trees_dbh.txt")
    cat("Reading in stemCloud from", clustList.name, "\n")
    if (!file.exists(clustList.name)) {
      cat("\nThere was no trees_dbh.txt found!\n")
      cat("Ending routine... Make sure there is a tree cluster file.\n")
      return()
    }
    clustList <- read.csv(clustList.name, sep = "\t")
    cat("-> Found", length(clustList[, 1]), "seed trees! ")

    # colnames(clustList)[colnames(clustList)=="cluster"] <- "id"

    idlist <- unique(clustList$id)
    cat("(Equ.", length(idlist), "unique tree ids?)\n")
    if (length(idlist) != length(clustList[, 1])) {
      warning("Careful! There are duplicates in the tree id-s.\n")
      cat("Assigning new tree ids...")
      clustList$oldId <- clustList$id
      clustList <- clustList[order(clustList$dbh, decreasing = T), ]
      clustList$id <- c(1:length(clustList$id))
      cat("Also dropping the z-coordinates and randomCol...")
      clustList$oldZ <- clustList$z
      clustList$oldRc <- clustList$randomCol
      clustList <- clustList[, -which(colnames(clustList) == "z")]
      clustList <- clustList[, -which(colnames(clustList) == "randomCol")]
      idlist <- unique(clustList$id)
      cat("Now we have", length(idlist), "unique tree ids!\n")
    }




    totalCloud <- add_lasattribute_manual(totalCloud, 0, "StemID", "Single Stem ID", type = "short")

    updateTreeList <- FALSE # writing out a new tree list with all the details z, dbh, randomcol

    allCols <- colnames(clustList)


    # FIXING Z MISSINGNOS
    if (!is.element("z", allCols) | sum(is.na(clustList$z)) > 1 | sum(clustList$z == 0.0) >= 10) {
      updateTreeList <- TRUE
      cat("\nAttaching missing z-values from ground model... ")
      dtmFile <- paste0(dirPath, groundPath, fileFinder, "_ground_min.grd")
      if (!file.exists(dtmFile)) {
        dtmFile <- paste0(dirPath, groundPath, fileFinder, "_ground_clip.grd")
      }
      tryCatch(
        {
          # read in raster file
          dtm_z <- raster(dtmFile)
        },
        error = function(error_condition) {
          cat("Error in reading the dtm-model, next loop!")
          next
        }
      )
      clustList$z <- round(extract(dtm_z, SpatialPoints(data.frame(x = clustList$x, y = clustList$y))) + 1.3, 3)
      # extracting the z-value from a dtm, was hard work to find that out... so easy
      # extract(dtm_c, SpatialPoints(data.frame(x=15, y=-10)))
      cat("done!\n")
    }

    # FIXING RANDOMCOL MISSINGNOS
    # if(is.element("randomCol", allCols)){
    #   clustList$randomCol[is.na(clustList$randomCol)] <- ""
    # }
    if (is.element("species", allCols)) {
      clustList$species[is.na(clustList$species)] <- ""
    }

    if (class(clustList$randomCol) == "character") {
      cat("Moving the comments from randomCol to a separate column!\n")
      clustList$comment <- clustList$randomCol
      clustList$randomCol <- NULL
      allCols <- colnames(clustList)
    }
    if (!is.element("randomCol", allCols) | sum(clustList$randomCol == 0) >= 2) {
      updateTreeList <- TRUE
      cat("We need new random numbers, as there were none in the stem list!\n")
      rnlist <- sample(1:65535, length(idlist))
      clustList$randomCol <- rnlist
    }


    findDBH <- FALSE
    # FIXING DBH MISSINGNOS
    ### NEW DBH MEASUREMENT FROM COORDINATES
    if (!is.element("dbh", allCols) || sum(is.na(clustList$dbh)) > 0) {
      findDBH <- TRUE
      cat("New DBH measurements for seed point cylinders will be done!")
      if (!is.null(clustList$dbh)) {
        cat("Retaining old dbh measurements.\n")
        clustList$dbh_old <- clustList$dbh
      }
      clustList$dbh <- 80
    }
    if (sum(clustList$dbh == 80) == length(clustList$dbh)) {
      cat("We have a old corrupt DBH measurement, so we are doing new ones!")
      findDBH <- TRUE
    }


    # sorting columns to x y z id randomCol dbh
    col_order <- c("x", "y", "z", "id", "randomCol", "dbh")
    copyList <- clustList
    clustList <- clustList[, c(col_order)]
    clustList <- cbind(clustList, copyList)
    clustList <- clustList[, !duplicated(colnames(clustList))]

    if (updateTreeList) {
      oldFile <- paste0("trees_dbh_", format(Sys.time(), "%Y%m%d_%H%M"), "_input.txt")

      if (file.rename(from = clustList.name, to = paste0(dbhPath, oldFile))) {
        cat("File", oldFile, "succesfully created.\n")
      } else {
        cat("Something went wrong with renaming sophisticated trees_dbh.txt (file was not present)...\n")
        return()
      }
      write.table(clustList,
        file = clustList.name,
        row.names = FALSE, sep = "\t"
      )
    }
    # 22-05-31: Cutting away the outer areas
    # plot(clustList$y ~ clustList$x, asp = 1)
    # abline(v = c(totalCloud@header@PHB$`Max X`, totalCloud@header@PHB$`Min X`), col = "red")
    # abline(h = c(totalCloud@header@PHB$`Max Y`, totalCloud@header@PHB$`Min Y`), col = "red")
    clustList <- clustList[!is.na(clustList$z), ]
    clustList <- clustList[clustList$x >= totalCloud@header@PHB$`Min X` & clustList$x <= totalCloud@header@PHB$`Max X`, ]
    clustList <- clustList[clustList$y >= totalCloud@header@PHB$`Min Y` & clustList$y <= totalCloud@header@PHB$`Max Y`, ]

    cat("Reduced the clustList roughly from", length(idlist), "to", length(clustList$id), "trees in", fileFinder, "\n")
    idlist <- unique(clustList$id)




    # 04-29: Random colors are now implemented in g1_ at the clustering process

    # 04-30 loading random colors from clustList file
    tryCatch(
      {
        rnlist <- clustList$randomCol
      },
      error = function(e) {
        cat("Something went wrong with reading random stem numbers...\n")
        return("PROBLEM IN RANDOM NUMBERS!")
        # cat("It is not good, but we are creating new random colors now...\n")
        # rnlist <- sample(1:65535, length(idlist))
      }
    )
    totalCloud <- add_lasattribute_manual(totalCloud, 0, "randomCol", "Random Stem ID", type = "short")

    if (mergedSet) totalCloud <- add_lasattribute(totalCloud, 0, "size", "Cluster size")
    # lidR::plot(stemCloud)
    cat("\n")





    if (useTreeFile != "") {
      if (file.exists(paste0(dbhPath, "seedLAS_cylinders.las"))) {
        oldFile <- paste0(dbhPath, "seedLAS_cylinders_before", format(Sys.time(), "%Y%m%d_%H%M"), ".las")
        file.rename(
          from = paste0(dbhPath, "seedLAS_cylinders.las"),
          to = oldFile
        )
      }
    }
    if (file.exists(paste0(dbhPath, "seedLAS_cylinders.las"))) {
      cat("Reading old seedLAS_cylinders.las file... ")
      seedLAS <- readLAS(paste0(dbhPath, "seedLAS_cylinders.las"), select = selector)
      seedClusts <- length(unique(seedLAS@data$StemID))
      cat("- ", seedClusts, "stemIDs found!\n")
    } else {
      # cut out a circle from point cloud
      # frame.up <- 0.3 #m
      # frame.down <- 0.3 #m
      # frame.rad <- 1.6 #times DBH


      cat("Cutting a cylinder of each seed tree:\n")
      # for every tree! put into seed set
      i <- 1
      seedLAS <- filter_poi(totalCloud, StemID != 0) # getting every detected stempoint as seed

      if (tileClipping == 0) {
        cat(
          "Efficiency improvement: Cutting slice from", min(clustList$z, na.rm = TRUE), "to",
          max(clustList$z, na.rm = TRUE), "m.\n"
        )
        cat("Reducing total", thMk(totalCloud@header@PHB$`Number of point records`), "points to ")
        reducedTotalCloud <- filter_poi(
          totalCloud, Z > min(clustList$z, na.rm = TRUE) - frame.down,
          Z < max(clustList$z, na.rm = TRUE) + frame.up
        )
        cat("cut", thMk(reducedTotalCloud@header@PHB$`Number of point records`), "points.\n")
      }

      start <- Sys.time()
      dropoutList <- clustList[0, ]

      if (findDBH) {
        circlePath <- paste0(crownPath, "circleFits/")
        if (!dir.exists(circlePath)) dir.create(circlePath)
      }
      clustListSafe <- clustList

      grabDBH <- function() {
        cat("Grab DBH function is not implemented! Refer to additional script allgo_addon_measure_dbh.R (useful)")
      }

      ## TILE GENERATION
      if (tileClipping > 0) {
        tileNumber <- tileClipping # 5x5 = 25
        tileBuffer <- 1

        xMax <- max(clustList$x) + 0.1
        xMin <- min(clustList$x) - 0.1
        yMax <- max(clustList$y) + 0.1
        yMin <- min(clustList$y) - 0.1

        tileSize.x <- xMax - xMin
        tileSize.y <- yMax - yMin

        gridX <- tileSize.x / tileNumber
        gridY <- tileSize.y / tileNumber



        cat("Spanning net of", tileNumber * tileNumber, "tiles x", tileSize.x, "and y", tileSize.y, "m (buffer is", tileBuffer, "m).\n")
        cat(" -> Single tile has x", gridX, "and y", gridY, "m.\n")

        xLs <- seq(from = xMin, to = xMax, by = gridX)
        yLs <- seq(from = yMin, to = yMax, by = gridY)

        alltexts <- merge(xLs, yLs)
        gridTileName <- function(x, y, noLetters = FALSE) {
          ydif <- yMax - y
          xdif <- x - xMin

          number <- trunc(xdif / gridX)
          row <- trunc(ydif / gridY)
          outName <- (number + tileNumber * row)
          return(outName)
        }


        alltexts$tileName <- (gridTileName(alltexts$x + 0.01, alltexts$y - 0.01, noLetters = TRUE))


        # clip(x1 = XU, x2 = XO, y1 = YU, y2 = YO)
        # abline(h = c(YO, YU), v = c(XO, XU), col = "blue")


        pdf(file = paste0(dbhPath, "tiles", tileClipping, "x", tileClipping, ".pdf"), width = 7, height = 7)
        tryCatch(
          {
            cexFactor <- 5 / mean(gridX + gridY)

            plot(clustList$x, clustList$y, type = "n", xlim = c(xMin, xMax), ylim = c(yMin, yMax), asp = 1)
            clip(x1 = xMin, x2 = xMax, y1 = yMin, y2 = yMax)
            abline(h = yLs, v = xLs, col = "grey")
            text(
              x = alltexts$x + gridX / 2,
              y = alltexts$y - gridY / 2,
              labels = alltexts$tileName, cex = 1.2
            )
            text(
              x = clustList$x + 0.5 + cexFactor,
              y = clustList$y - 0.5 - cexFactor,
              labels = clustList$id, cex = cexFactor * 2
            )
            points(clustList$x, clustList$y, cex = cexFactor)

            dev.off()
          },
          error = function(e) {
            dev.off()
            cat("Error when writing the .pdf, closing device.\n")
          }
        )


        clustList$tileName <- strtoi(gridTileName(clustList$x, clustList$y, noLetters = TRUE), base = 10)
        # hist(clustList$tileName)

        tiles <- unique(clustList$tileName)
        tiles <- tiles[order(tiles)]

        cat("\nProcessing seeds (frame.rad =", frame.rad, "x):\n")
        for (a in 1:length(tiles)) {
          if (a %% tileClipping == 0) {
            rm(tilesLAS)
            gc()
          }
          cat("-- tile", tiles[a])
          # sta <- Sys.time()

          clustList.sub <- clustList[clustList$tileName == tiles[a], ]

          xu <- alltexts$x[alltexts$tileName == tiles[a]]
          xo <- xu + gridX
          yo <- alltexts$y[alltexts$tileName == tiles[a]]
          yu <- yo - gridY

          tilesLAS <- filter_poi(
            totalCloud, X > xu - tileBuffer, X < xo + tileBuffer,
            Y > yu - tileBuffer, Y < yo + tileBuffer
          )
          cat("", tilesLAS@header@PHB$`Number of point records`, "pts")


          for (b in 1:length(clustList.sub$id)) {
            i <- b
            if (b == 1) {
              cat(" - ")
            } else {
              cat(", ")
            }
            cat(clustList.sub$id[b])
            # sta <- Sys.time()
            nowClust <- clustList.sub[b, ]
            # print(nowClust)
            nowZ <- nowClust$z
            # cat("b")

            cutRadius <- nowClust$dbh
            if (cutRadius <= 5) {
              cutRadius <- 5
            }
            clustCloud <- clip_circle(tilesLAS,
              xcenter = nowClust$x, ycenter = nowClust$y,
              radius = cutRadius / 200 * frame.rad
            )
            # cat("c", clustCloud@header@PHB$`Number of point records`)
            # cat("Reduced seed to", thMk(clustCloud@header@PHB$`Number of point records`), "points.\n")
            clustCloud <- filter_poi(clustCloud, Z > nowClust$z - frame.down, Z < nowClust$z + frame.up)

            if (clustCloud@header@PHB$`Number of point records` == 0) {
              cat(" DISCARD: no seed points at DBH found (!)")
              dropoutList[length(dropoutList[, 1]) + 1, ] <- nowClust
              next
            }

            if (findDBH & nowClust$id >= 9000) {
              # grabDBH()
            }


            # cat("d", clustCloud@header@PHB$`Number of point records`)
            clustCloud@data$StemID <- nowClust$id

            # cat("e", clustCloud@header@PHB$`Number of point records`)
            clustCloud@data$randomCol <- nowClust$randomCol
            # cat(clustCloud@header@PHB$`Number of point records`, "points. - -\n")
            seedLAS@data <- rbind(clustCloud@data, seedLAS@data)
            # sto <- Sys.time()
            # print.difftime(sto-sta)
          }
          cat("\n\n")
        }
      } else {
        # clip every single tree individually to seed unit
        for (i in 1:length(clustList[, 1])) {
          cat(" - seed", i)
          # sta <- Sys.time()
          nowClust <- clustList[i, ]
          # print(nowClust)
          nowZ <- nowClust$z
          # cat("b")

          clustCloud <- clip_circle(reducedTotalCloud,
            xcenter = nowClust$x, ycenter = nowClust$y,
            radius = nowClust$dbh / 200 * frame.rad
          )
          # cat("c", clustCloud@header@PHB$`Number of point records`)
          # cat("Reduced seed to", thMk(clustCloud@header@PHB$`Number of point records`), "points.\n")
          clustCloud <- filter_poi(clustCloud, Z > nowClust$z - frame.down, Z < nowClust$z + frame.up)

          if (clustCloud@header@PHB$`Number of point records` == 0) {
            cat(" DISCARD: no seed points at DBH found (!)")
            dropoutList[length(dropoutList[, 1]) + 1, ] <- nowClust
            next
          }

          if (findDBH & nowClust$id >= 9000) {
            # grabDBH()
          }


          # cat("d", clustCloud@header@PHB$`Number of point records`)
          clustCloud@data$StemID <- nowClust$id

          # cat("e", clustCloud@header@PHB$`Number of point records`)
          clustCloud@data$randomCol <- nowClust$randomCol
          # cat(clustCloud@header@PHB$`Number of point records`, "points. - -\n")
          seedLAS@data <- rbind(clustCloud@data, seedLAS@data)
          # sto <- Sys.time()
          # print.difftime(sto-sta)
        }
      }






      if (length(dropoutList[, 1]) != 0) {
        cat("\nThese seeds were deleted as false positives, no points present at DBH region:\n")
        print(dropoutList)
        clustList <- clustList[!is.element(clustList$id, dropoutList$id), ]
        write.table(dropoutList,
          file = paste0(dbhPath, "trees_dropped.txt"),
          row.names = FALSE, sep = "\t"
        )
      } else {
        cat("\n")
      }

      # always do that, because we do want to have a
      # honest split set of not-included trees in the current folder
      {
        cat("Writing new trees_dbh.txt without all border trees DBHs... \n")
        oldFile <- paste0("trees_dbh_all.txt")

        if (file.rename(from = clustList.name, to = paste0(dbhPath, oldFile))) {
          cat("File", oldFile, "succesfully renamed, and new trees_dbh.txt created.\n")
        } else {
          cat("Something went wrong with renaming all trees_dbh.txt (file was not present)...\n")
          return()
        }
        write.table(clustList,
          file = clustList.name,
          row.names = FALSE, sep = "\t"
        )
      }


      cat("Generating seedLAS_cylinders.las... \n")
      writeLAS(seedLAS, paste0(dbhPath, "seedLAS_cylinders.las"))
      seedClusts <- length(unique(seedLAS@data$StemID))
      cat("  Check: There are", seedClusts, "seeds and", length(unique(clustList$id)), "trees in the list")


      stop <- Sys.time()
      print.difftime(round(stop - start, 1))
    }






    tallyTrees <- sum(clustList$dbh > 5)
    cat("There are", tallyTrees, "trees larger than 5 cm\n")
    if (tallyTrees < 0.2 * length(idlist)) {
      cat("Only small trees in our set! Considering all trees for limitStems.\n")
      tallyTrees <- length(idlist)
    }
    minimumNowStems <- limitStems / 100 * tallyTrees
    cat("We need at least", minimumNowStems, "seed stems (", limitStems, "%) per round or we will increase.\n")








    if (!silent) lidR::plot(seedLAS)

    # CUTTING INTO WINDOW
    if (sum(cutWindow == c(-1000, -1000, 2000)) != 3) {
      cat("\nPre-cutting totalCloud:\n")
      cat("a) Reducing seedLAS from", thMk(seedLAS@header@PHB$`Number of point records`), "points ")
      seedLAS <- filter_poi(seedLAS, X > XL, X < XL + width, Y > YL, Y < YL + width)
      cat("to", thMk(seedLAS@header@PHB$`Number of point records`), "points.\n")

      cat("b) Reducing totalCloud from", thMk(totalCloud@header@PHB$`Number of point records`), "points ")
      totalCloud <- filter_poi(totalCloud, X > XL, X < XL + width, Y > YL, Y < YL + width)
      # frameLAS <- filter_poi(stemCloud, X > min(xrange), X < max(xrange), Y > min(yrange), Y < max(yrange))
      cat("to", thMk(totalCloud@header@PHB$`Number of point records`), "points.\n")
    }


    # VOXELISATION OF INPUT POINTS ####
    if (voxelSize != 0) {
      cat("Reduction into VOXEL with a cube length of", voxelSize, "cm... \n")

      seedLAS <- tlsSample(seedLAS, smp.voxelize(voxelSize / 100))
      # TODO: Change that into normal voxels, but then all the rest of the script goes down...
      cat(" VOXEL POINTS in total ")
      cat(thMk(seedLAS@header@PHB$`Number of point records`), "seed points and ")

      allocationCloud <- totalCloud

      groundCloud <- voxelize_points(groundCloud, voxelSize / 100)
      groundCloud@data$Classification <- 2L
      cat(thMk(groundCloud@header@PHB$`Number of point records`), "ground points and ")

      vegCloud <- voxelize_points(vegCloud, voxelSize / 100)
      vegCloud@data$Classification <- 0L
      cat(thMk(vegCloud@header@PHB$`Number of point records`), "vegetation points remain.\n")
      cat("Joining them with proper Classification (2=ground 0=veg)... ")
      totalCloud <- rbind(groundCloud, vegCloud)
      rm(groundCloud, vegCloud)
      gc()
      cat("done!\n")
    }



    # length(stemCloud@data$X)
    cat("Adding seed points... ")
    stemCloud <- totalCloud
    stemCloud@data <- rbind(seedLAS@data, stemCloud@data, fill = TRUE)
    stemCloud@data$StemID[is.na(stemCloud@data$StemID)] <- 0
    stemCloud@data$Classification[is.na(stemCloud@data$Classification)] <- 2L

    if (is.element("gpstime", colnames(stemCloud@data))) {
      stemCloud@data$gpstime[is.na(stemCloud@data$gpstime)] <- 555
    }

    if (is.element("R", colnames(stemCloud@data))) {
      cat("(dropping the RGB information")
      removes <- c(
        which(colnames(stemCloud@data) == "R"),
        which(colnames(stemCloud@data) == "G"),
        which(colnames(stemCloud@data) == "B")
      )
      subData <- subset(stemCloud@data, select = -c(removes))
      stemCloud <- LAS(subData)
      cat(") ")
    }



    pointsBefore <- sum(is.na(stemCloud@data$StemID))
    stemCloud <- filter_duplicates(stemCloud)
    pointsAfter <- sum(is.na(stemCloud@data$StemID))

    if (pointsBefore != pointsAfter) {
      cat("PROBLEM: WE LOST SOME SEED POINTS BY DUPLICATE FILTERING!!! \n\n\n")
    }
    # cat(thMk(stemCloud@header@PHB$`Number of point records`), "points in total.\n")

    # writeLAS(stemCloud, "stems.las")
    stemCloud <- add_lasattribute_manual(stemCloud, 999, "runJay", "number of run added to total", type = "short")
    rm(totalCloud)
  }

  ##### PREPARATION #####
  # cutting out sample rectangle
  if (sum(cutWindow == c(-1000, -1000, 2000)) == 3) {
    cat("No cutting because of no cutWindow setting!\n")
    frameLAS <- filter_poi(stemCloud, Intensity > -1000000000)
  } else {
    cat("Cutting to the dimensions... ")
    frameLAS <- filter_poi(stemCloud, X > XL, X < XL + width, Y > YL, Y < YL + width)
    # frameLAS <- filter_poi(stemCloud, X > min(xrange), X < max(xrange), Y > min(yrange), Y < max(yrange))
    cat(thMk(frameLAS@header@PHB$`Number of point records`), "points remain.\n")
  }
  rm(stemCloud)
  ### HUOM: Originally I wanted to voxelize the final input cloud right before region growing
  # in fact, the mixture of seed points and blank points results in a loss of classification (and all the other values)
  # take care about it in the end also when outputting (one last big round of neighbour allocations
  # (all voxel points color all raw input points))
  # now voxelisation is done at the creation of the seed las set and right before merging both
  if (zScale != 1) {
    frameLAS@data$Z <- frameLAS@data$Z / zScale
  }
  # lidR::plot(frameLAS, color = "StemID")
  pointNumber <- frameLAS@header@PHB$`Number of point records`
  cat("Framed file containing:", thMk(pointNumber), "points.\n")


  if (doReferencedOnly) {
    cat("We reset out the", thMk(sum(frameLAS@data$StemID > 20000)), "bushy seeds over 20000 ID to blanks...\n")
    # here we assign the > 20.000er IDs to blanks again, we dont make trees out of those crappy hits
    frameLAS@data$StemID[frameLAS@data$StemID > 20000] <- 0
  }

  seeds <- filter_poi(frameLAS, StemID != 0) # getting every detected stempoint as seed
  seeds@data$runJay <- 0
  # head(seeds@data[, c(1:3, 20)])
  cat(
    "We remain starting with", thMk(seeds@header@PHB$`Number of point records`), "seed points (=",
    round(seeds@header@PHB$`Number of point records` / pointNumber * 100, 1), "%)\n"
  )


  ### EFFICIENCY SET HERE, LASCATALOG MORE EFFICIENT WITH NEIGHBORS!!! 04-30
  seedSet <- seeds@data # for faster handling
  if (!silent) lidR::plot(seeds, color = "StemID")


  if (mergedSet) { # if merged set, prefer the smaller ones by sorting and duplicate removing
    seedSet <- seedSet[order(seedSet$size, decreasing = FALSE), ]
    # head(seedSet[, c(1:3, 20, 21)])
    seeds@data <- seedSet
    seeds <- filter_duplicates(seeds)
  }

  outLAS <<- seeds # gathers all assigned points

  blankLAS <<- filter_poi(frameLAS, StemID == 0) # gathers all blank points
  blankSetLength <- length(blankLAS@data$X)

  if (!silent) {
    # writing out all the input file without segmentation
    if (zScale != 1) {
      frameLAS@data$Z <- frameLAS@data$Z * zScale
    }
    writeLAS(frameLAS, paste0(crownPath, fileFinder, "_input_all", locationStr, ".laz"))
  }
  rm(frameLAS)
  gc()

  distanceCounter_cm <- 0
  # 22-06-08: distanceCounter can be set at the top
  # distanceCounter_cm_limit <- 200
  groundFound <- FALSE

  useableDistance <- 0.7
  justUp <- FALSE # if now incrementing distance, next time take all in!

  #### watch out, different concept UP != DOWN ;)
  # justDown <- TRUE #this variable remembers if instantly went down the searching distance,
  # if that happened and next round we need to go up again, no taking of all seeds again.
  # NOT IMPLEMENTED: leaving out some seeds as they are not growing (stagnating)


  # old settings normal
  start.from <- 0.15
  start.to <- 0.01
  start.by <- 0.0025
  # try a new
  start.from <- 0.10
  start.to <- 0.01
  start.by <- 0.005
  # try a new 2
  start.from <- 0.01
  start.to <- 0.03
  start.by <- 0.005

  if (voxelSize != 0) {
    start.from <- voxelSize / 100 + 0.0001
    start.to <- (voxelSize) / 100 + 0.0001 # stays the same, voxelsize steplength of 1 for voxelised growing

    if (diagonals) {
      cat("We use initial cube diagonals, so we increase the search distance by sqrt(2).\n")
      start.from <- round(voxelSize / 100 * 1.6, 2) + 0.0001
      start.to <- round(voxelSize / 100 * 1.6, 2) + 0.0001
    }
    start.by <- 0.005
    incrementDistance <- voxelSize / 100
    cat("Due to voxelisation, the increment is set to one voxel or", voxelSize, "cm.\n")


    start.times <- ceiling(distanceCounter_cm_limit / voxelSize / zScale)
    last.quart <- floor(start.times / 4) # the last fourth * sqrt 2 to get the diagonals as well (roundish footprint)

    startRow <- c(
      rep(round(voxelSize / 100, 2) + 0.0001, start.times - last.quart),
      rep(round(voxelSize / 100 * 1.6, 2) + 0.0001, last.quart)
    ) # times sqrt(2) for circleish shape on the ground
    cat(paste0(start.times, " Rounds of the start row: ", start.times - last.quart, "x ", round(voxelSize, 0), " cm (simple growth) and ", last.quart, "x ", round(voxelSize * 1.6, 0), " cm (diagonal gr.).\n"))
  } else {
    start_cm <- 2
    start.times <- distanceCounter_cm_limit / start_cm / zScale
    startRow <- rep(start_cm / 100, start.times)
    cat("We repeat", start.times, "times the start row at", start_cm / 100, "m\n")
  }


  reachableMaxHeight <- sum(startRow) * zScale * useableDistance
  cat(
    "Adding the fast climbing", reachableMaxHeight,
    "m to the seed-maximum height of", round(max(seeds@data$Z) * zScale, 1),
    "m we should be at approx.", round(max(seeds@data$Z) * zScale + reachableMaxHeight, 1),
    "m in round", length(startRow), "\n"
  )

  cat("We are doing", totalRuns, "rounds in total, maximum", durationMins, "minutes.\n")





  #### START CROWN FEELING ####
  useAllBlanks <- FALSE # if 99 % of blanks are already in the set, we dont need to split them (time saving)
  moreZ <- 0.11 # initial value for z-cutting in blank set (11 cm higher than highest point, rest is cut off
  # this value is later set to an adaptive value according to searching distance)


  #### Crowns: loops ####
  prgsFrame <- data.frame(
    "round" = 0L, "dist" = 0, "cumDist" = 0, "maxh" = 0, "pctAdded" = 0, "treesGrow" = 0,
    "blanks" = 0, "blanks.used" = 0, "seeds" = 0,
    "pctBlanks" = 0, "pctSeeds" = 0, "pctPuzzle" = 0,
    "timeNB" = 0, "timeTotal" = 0
  )


  despr <- 0 # take in all again only if 10x (despr.max) no success
  despr.max <- 5
  decreas.count <- 0
  decreas.max <- 20 # if 20 times larger than adding minimum, decreasing search distance

  startRowFinished <- FALSE # if once too less points, we stop the decreasing run

  # search distance explanation - how is it calculated?
  # starting with 1 cm (special case running down from 15 cm in the first rounds 15 14 13 usw.)
  # then it is only increased by the number of incrementDistance (0.5cm),
  # if there are less than 0.5 % of total blanks that are assigned now
  # this variable is called "limitShare"
  # every time this happens, all seed points are taken in again (more distance - could find new points)
  # if the distance is larger than 0.05, so 5 cm, then not always all are taken in again,
  # only if tha increment has happened the 10th time (setting in despr. max)
  # additionally, every 20th time (decreas.max) the setting has worked out fine,
  # the searching distance is decreased again by incrementDistance
  # for every increase, the decreasing count is punished with an additional 10 rounds to keep stable


  # limitShare <- 0.005 # from function header
  # that variable says how much of all blank points must at least be assigned in next step
  # if less than half a percent seed points left, reset seeding to all
  # incrementDistance <- 0.005 # 5 mm more per increment
  # maximumDistance <- 0.50 # if nearest points are further than 50 cm, for is terminated early





  for (j in 1:totalRuns) {
    start <- Sys.time()
    cat("\n")


    if (j == 1) {
      searchDistance <- startRow[1] # 18cm start, first 15 rounds steadily decreased
      cat(fileFinder, ": Iteration ", sprintf("%03d", j),
        " on _", format(Sys.time(), "%a%d. %X"), "_",
        sep = ""
      )

      # RUN 1 - special
      maxZ <- max(seeds@data$Z) + (moreZ / zScale)
      cat("max height", round(maxZ * zScale, 2), "m, cumDist 0 cm\n")


      blankSetLength <- nrow(blankLAS@data)
      blanksUsed <- 100
      if (!useAllBlanks) {
        blankSet.sub <- filter_poi(blankLAS, Z < (maxZ))@data
        blankLAS <<- filter_poi(blankLAS, Z >= (maxZ))
        blanksUsed <- round(nrow(blankSet.sub) * 100 / (blankSetLength), 1)
        cat("blanks in race:", thMk(nrow(blankSet.sub)), "(", blanksUsed, "% of all blanks)\n")
      } else {
        blankSet.sub <- blankLAS@data
        cat("blanks in race:", thMk(blankSetLength), "(ALL)\n")
      }

      cat("and seeds in race:", thMk(nrow(seedSet)), "(", round(nrow(seedSet) * 100 / (pointNumber), 2), "% )\n")

      # cat("blanks:", thMk(blankSetLength), "(", round(blankSetLength*100 / (pointNumber), 1), "% )")
      # cat(", seeds:", thMk(length(seedSet$X)), "(", round(length(seedSet$X)*100 / (pointNumber), 2), "% )\n")
      # cat("blanks in race:", thMk(length(blankSet.sub$X)), "(", blanksUsed, "% of all blanks)\n")

      cat("nn2... ")
      # closest <- nn2(data.frame(seedSet[, 1:3]), query = data.frame(blankSet.sub[, 1:3]),
      #                k = howManyClosestPoints, searchtype = "standard")
      select_cols <- c("X", "Y", "Z")
      closest <- as.data.table(
        nn2(seedSet[, ..select_cols],
          query = blankSet.sub[, ..select_cols],
          k = howManyClosestPoints, searchtype = "radius", radius = searchDistance
        )
        #       k = howManyClosestPoints, searchtype = "standard", radius = searchDistance)
      )
      stop <- Sys.time()
      # cat("Closest neighbors found.\n")
      timeNB <- as.difftime(stop - start)
      start <- Sys.time()
    } else {
      # RUN 2 - 1000
      start <- Sys.time()

      # NEW RULE HERE 8 / 20, no caring for minimum number of stems, if half of maximum distance is already reached! or it will soon be over, expecially with voxel.
      if (length(nowStems) < minimumNowStems && searchDistance < maximumDistance / 2) {
        cat("Taking all points in because less than", minimumNowStems, "stems added in last round!\n\n")
        searchDistance <- searchDistance + incrementDistance
        despr <- 0
        gc() # garbage collection
        seedSet <- outLAS@data
      } else {
        # Setting distance for searching and seed set

        if (nowShare < limitShare * 100 && startRowFinished) {
          decreas.count <- 0
          searchDistance <- searchDistance + incrementDistance
          despr <- despr + 1
          # alernaternative with justdown: condition = (despr > despr.max || (searchDistance < 0.05 && !justDown))
          if ((justUp && nrow(blankLAS@data) / pointNumber > 0.05) || despr > despr.max || searchDistance < 0.05) { # this is kind of a border ...| 5cm |.... no jumps
            if (justUp && voxelSize < 1) {
              justUp <- FALSE # now we did the just up extra handling, dont want it a second time!
              searchDistance <- searchDistance + incrementDistance
              # incrementing again (Double effect)
            } else {
              justUp <- TRUE # only if not just taken all in again
            }
            ## bad punishment...
            despr <- 0
            gc() # garbage collection
            cat("Too less seed points, taking all in again...\n\n")
            seedSet <- outLAS@data
          } else {
            justUp <- TRUE
            seedSet <- seeds@data
            if (voxelSize > 1) {
              despr <- 0
              gc() # garbage collection
              cat("Undershooting limitshare in Voxel settings, taking all in again...\n\n")
              seedSet <- outLAS@data
            }
            # seedSet <- filter_poi(outLAS, runJay > j-3)@data # never tried that, only using 3 last times seeds
            # could errorize because some seeds are lost forever, need only sometimes last ones, sometimes all
          }

          if (searchDistance > maximumDistance) {
            cat("Terminating search - maximum distance exceeded in run", j, "\n\n")
            break
          }

          cat(fileFinder, ": Iteration ", sprintf("%03d", j),
            " on _", format(Sys.time(), "%a%d. %X"), "_     ",
            sep = ""
          )
        } else {
          cat(fileFinder, ": Iteration ", sprintf("%03d", j),
            " on _", format(Sys.time(), "%a%d. %X"), "_     ",
            sep = ""
          )


          if (j <= length(startRow)) {
            cat("STROW ")
          }
          justUp <- FALSE
          # justDown <- FALSE
          seedSet <- seeds@data
          if (nowShare > 1.5 * limitShare * 100) {
            decreas.count <- decreas.count + 1
          }

          # old condition: if(decreas.count >= decreas.max)
          if (nowShare > 1.5 * limitShare * 100 && decreas.count >= decreas.max) {
            # FINALLY DO THE DECREASING!
            decreas.count <- 0
            searchDistance <- searchDistance - incrementDistance
            cat("DROP ")
            # justDown <- TRUE # just decreased if we are instantly going up again
            # one more chance is from the beginning and if we are so stable that we could decrease distance search
          }
        }
      }

      if (voxelSize != 0) {
        if (searchDistance < voxelSize / 100) searchDistance <- voxelSize / 100 + 0.0001
      }


      moreZ <- useableDistance * searchDistance * zScale # adaptive cutting off only for next run
      # 70 % before ensures an equal height progress of all trees (highest tree cannot grow away)
      # if you want a single tree specific height progress (individual growing) use following line:
      # moreZ <- searchDistance * zScale # every tree can grow to maximum extend maxdist in zScale
      # warning: if enabled, highest trees could outgrow smaller ones if search distance is getting large
      # additionally slower for large files: the big trees require the whole vegetation cloud (only max Z of total cloud is cut)
      # but: if there are steep slopes it makes sense to higher the useable distance, or the highest trees might be stuck small



      if (distanceCounter_cm * zScale > distanceCounter_cm_limit &&
        !startRowFinished) {
        startRowFinished <- TRUE
        if (voxelSize != 0) {
          cat("RESET ")
          searchDistance <- startRow[1]
        }
      }

      # 200 tries to assign every point to a stem
      if (j < length(startRow) && !startRowFinished) {
        searchDistance <- startRow[j]
        moreZ <- useableDistance * searchDistance * zScale # adaptive cutting off only for next run
      }





      cat(
        "dist =",
        round(searchDistance * 100), "cm\n"
      )



      # lidR::plot(seeds, color = "StemID")
      # lidR::plot(blankLAS, color = "StemID")
      # lidR::plot(rbind(blankLAS, blanks), color = "StemID")

      if (blanksUsed > 80) useAllBlanks <- TRUE
      maxZ <- max(seeds@data$Z) + (moreZ / zScale)
      cat("max height", round(maxZ * zScale, 2), "m, cumDist", distanceCounter_cm, "cm, ")
      cat("done:", thMk(outLAS@header@PHB$`Number of point records`), "(", round(outLAS@header@PHB$`Number of point records` * 100 / (pointNumber), 2), "% )\n")

      blankSetLength <- nrow(blankLAS@data)
      if (!useAllBlanks) {
        blankSet.sub <- filter_poi(blankLAS, Z < (maxZ))@data
        blankLAS <<- filter_poi(blankLAS, Z >= (maxZ))
        blanksUsed <- round(nrow(blankSet.sub) * 100 / (blankSetLength), 1)

        if (length(blankSet.sub) == 0) {
          useAllBlanks <- TRUE
          blankSet.sub <- blankLAS@data
          cat("blanks in race:", thMk(blankSetLength), "(ALL)\n")
        } else {
          cat("blanks in race:", thMk(nrow(blankSet.sub)), "(", blanksUsed, "% of all blanks)\n")
        }
      } else {
        blanksUsed <- 100
        blankSet.sub <- blankLAS@data
        cat("blanks in race:", thMk(blankSetLength), "(ALL)\n")
      }
      cat("and seeds in race:", thMk(nrow(seedSet)), "(", round(nrow(seedSet) * 100 / (pointNumber), 2), "% )\n")



      # blanksUsed <- round(nrow(blankSet.sub)*100 / (blankSetLength), 1)
      # if(!useAllBlanks) cat("using only:", blanksUsed, "% of left unassigned points -> cut to ", thMk(length(blankSet.sub$X)), ".\n")
      if (blankSetLength == 0) {
        cat("Terminating searching - all blanks allocated!\n")
        break
      }



      cat(" nn2... ")
      # closest <- nn2(data.frame(seedSet[, 1:3]), query = data.frame(blankSet.sub[, 1:3]),
      #                k = howManyClosestPoints, searchtype = "standard")
      select_cols <- c("X", "Y", "Z")
      closest <- as.data.table(
        nn2(seedSet[, ..select_cols],
          query = blankSet.sub[, ..select_cols],
          k = howManyClosestPoints, searchtype = "radius", radius = searchDistance
        )
        #         k = howManyClosestPoints, searchtype = "standard")
      )
      # to do: add the round number j to information gathering more efficiently
      cat("found! ")
      stop <- Sys.time()
      timeNB <- as.difftime(stop - start)
      start <- Sys.time()
    }




    prgsFrame[j, ]$round <- j
    prgsFrame[j, ]$maxh <- round(maxZ * zScale, 2)
    prgsFrame[j, ]$blanks <- blankSetLength
    prgsFrame[j, ]$blanks.used <- blanksUsed
    prgsFrame[j, ]$pctBlanks <- round(blankSetLength * 100 / (pointNumber), 1)
    prgsFrame[j, ]$seeds <- nrow(seedSet)
    prgsFrame[j, ]$pctSeeds <- round(nrow(seedSet) * 100 / (pointNumber), 2)


    # add an point identifier for all Million seeds
    seedSet <- seedSet[, nn.idx := c(1:nrow(seedSet))]

    # assign StemID to every point
    closest <- merge(closest, seedSet[, c("StemID", "nn.idx")],
      all.x = TRUE, sort = FALSE
    )

    # remove all points that are too far away
    closest[nn.dists > searchDistance, c("StemID")] <- 0

    # keep_close_select <- closest$nn.dists <= (searchDistance/100)
    deletedPoints <- closest[StemID == 0, .N]
    assignedPoints <- (nrow(blankSet.sub) - deletedPoints)
    nowShare <- round(assignedPoints / blankSetLength * 100, 2)

    cat(
      "Assigning", thMk(assignedPoints), "pts of total",
      thMk(blankSetLength),
      "(nowShare =", nowShare, "%)\n"
    )



    distanceCounter_cm <- distanceCounter_cm + round(searchDistance * 100, 1)
    prgsFrame[j, ]$dist <- searchDistance
    prgsFrame[j, ]$cumDist <- distanceCounter_cm
    prgsFrame[j, ]$pctAdded <- nowShare
    prgsFrame[j, ]$pctPuzzle <- round(100 * assignedPoints / pointNumber, 3)
    # these numbers were used for troubleshooting with packaging and crown feeling problems

    # if(1 == 2){
    #
    #   # cat("1 -- ")
    #   stemidL <- list(StemID = rep(0L, length(blankSet.sub$X)))
    #
    #   # cat("2 -- ")
    #   for(i in 1:limit){
    #     stemidL$StemID[closest.sub$nn.point[i]] <- (seedSet$StemID[closest.sub$nn.idx[i]])
    #   }
    #
    #   # cat("3 -- ")
    #   blankSet.sub$StemID <- data.frame(stemidL)[, 1]
    #
    # }

    blankSet.sub$StemID <- closest$StemID




    # informations, about which seed tree entities have been assigned this round
    nowStems <- unique(blankSet.sub[StemID != 0, "StemID"])$StemID
    # nowStems <- nowStems[nowStems!=0] # 0 refers to the blanks, is no tree at all
    nowStems <- nowStems[order(nowStems)]

    cat(" - Our newly assigned blanks hold", length(nowStems), "tree entities ranging from", range(nowStems), "\n")
    if (length(nowStems) < 50) {
      cat(nowStems)
      cat("\n")
    }

    # cat("Merge outLAS... ")
    seeds@data <- blankSet.sub[StemID != 0]
    # cat("3 -- ")
    seeds@data$runJay <- j
    # cat("ok\n")

    blankSet.sub <- blankSet.sub[StemID == 0]
    # print((blankSet.sub$StemID))
    # cat(length(blankLAS@data$X), " x ", length(blankLAS@data), " und 0er len ", sum(blankSet.sub$StemID==0), " ")
    # cat(ncol(blankSet.sub[blankSet.sub$StemID==0]), " x ", nrow(blankSet.sub[blankSet.sub$StemID==0]))
    if (!useAllBlanks) {
      blankSet.sub <- rbind(blankSet.sub, blankLAS@data)
    }
    blankLAS@data <<- blankSet.sub
    # cat("2 -- ")


    # remove nn.idx column again
    outLAS@data <<- outLAS@data[, -c("nn.idx")]
    outLAS <<- rbind(outLAS, seeds)

    if (!silent) {
      if (j %% plotEvery == 1 || j < 3) {
        if (zScale != 1) {
          outLAS@data$Z <<- outLAS@data$Z * zScale
        }
        lidR::plot(outLAS, color = "StemID")
        if (zScale != 1) {
          outLAS@data$Z <<- outLAS@data$Z / zScale
        }
      }
    }

    stop <- Sys.time()
    timeAS <- as.difftime(stop - start)
    cat(
      "Elapsed time: NB", round(timeNB, 1), units(timeNB),
      "- AS", round(timeAS, 1), units(timeAS),
      "- TOTAL", round(timeAS + timeNB, 1), units(as.difftime(timeAS + timeNB)), "\n"
    )

    prgsFrame[j, ]$timeNB <- round(as.numeric(timeNB, units = "secs"), 1)
    prgsFrame[j, ]$timeTotal <- round(as.numeric(timeNB + timeAS, units = "secs"), 1)
    prgsFrame[j, ]$treesGrow <- length(nowStems)



    if (distanceCounter_cm * zScale > distanceCounter_cm_limit && !groundFound) {
      cat("\n")
      cat("We reached the z-limit (", distanceCounter_cm_limit, "cm), so we should be at the ground.\n")
      grLAS <- outLAS
      try(grLAS@data$Z <- grLAS@data$Z * zScale)
      grLAS <- add_lasattribute_manual(grLAS, grLAS@data$StemID, "StemID", "Single Stem ID", type = "short")
      grLAS <- add_lasattribute_manual(grLAS, 0, "randomCol", "Random Stem ID", type = "short")
      try(grLAS@data$randomCol <- rnlist[match(grLAS@data$StemID, idlist)])
      grLAS@data$randomCol[is.na(grLAS@data$randomCol)] <- 0
      writeLAS(grLAS, paste0(crownPath, fileFinder, "_grounded_stems", locationStr, ".las"))
      rm(grLAS)
      gc()
      reduceGroundPoints_num <- blankLAS@data[Classification == 2, .N]
      cat("Unassigned ground points:", thMk(reduceGroundPoints_num), "pts\n")
      pointNumber <- pointNumber - reduceGroundPoints_num
      cat("New reduced pointNumber:", thMk(pointNumber), "pts\n")
      blankLAS <<- filter_poi(blankLAS, Classification != 2)
      groundFound <- TRUE
      cat("Ground points are removed after round", j, "from region growing.\n\n")
    }






    if (as.numeric(Sys.time() - globalStart, units = "mins") > durationMins) {
      cat("Terminating loop because of time duration exceeded", durationMins, "minutes.\n\n")
      break()
    }
  }






  write.table(prgsFrame, paste0(crownPath, fileFinder, "_crowns_meta_", startTimeString, ".txt"),
    sep = "\t", row.names = FALSE
  )
  cat("Output meta table created in:\n   ", paste0(crownPath, fileFinder, "_crowns_meta_", startTimeString, ".txt"), "\n")


  cat("\nGenerating crown point cloud from outLAS:\n")
  outLAS <<- add_lasattribute_manual(outLAS, outLAS@data$StemID, "StemID", "Single Stem ID", type = "short")
  beforeLength <- outLAS@header@PHB$`Number of point records`
  cat("Voxel file contains in total", thMk(beforeLength), "pts.\n")

  # cat("After filtering...")
  # outLAS <<- filter_duplicates(outLAS)
  # afterLength <- outLAS@header@PHB$`Number of point records`
  # cat(" Now left are", thMk(afterLength), "points (removed", thMk(beforeLength - afterLength),
  #    "points =", round((beforeLength - afterLength) / beforeLength*100, 2), "%).\n")

  #### OUTPUTTING ####
  if (zScale != 1) {
    cat("Resetting z-Scale - ")
    if (blankSetLength != 0) {
      try(blankLAS@data$Z <<- blankLAS@data$Z * zScale)
    }
    try(outLAS@data$Z <<- outLAS@data$Z * zScale)
    cat("done.\n")
  }

  cat("Adding randomCol field... ")
  tryCatch(
    {
      if (!exists("rnlist")) {
        cat("We need new random numbers, as there were none in the stem list!\n")
        rnlist <- sample(1:65535, length(idlist))
      }

      outLAS <<- add_lasattribute_manual(outLAS, 0, "randomCol", "Random Stem ID", type = "short")
      outLAS@data$randomCol <<- rnlist[match(outLAS@data$StemID, idlist)]
      outLAS@data$randomCol[is.na(outLAS@data$randomCol)] <<- 0
    },
    error = function(e) cat("Something went wrong with random stem numbers...\n")
  )
  cat("done!\n")

  if (voxelSize != 0) {
    cat("\n")
    cat("(VoX-Back) Final allocation round on _", format(Sys.time(), "%a%d. %X"), "_:\n", sep = "")

    sta <- Sys.time()

    outLAS_name <<- paste0(crownPath, fileFinder, "_crowns", locationStr, ".laz")
    cat("Writing outLAS to /crowns_VOX.laz... ")
    writeLAS(outLAS, paste0(crownPath, fileFinder, "_crowns", locationStr, "_VOX.laz"))
    cat(thMk(outLAS@header@PHB$`Number of point records`), "pts ok.\n")

    # last allocation round of all voxelized points
    #
    # plot(outLAS) # output voxel points
    # allocationCloud #in fact totalcloud has all points

    allocationDistance <- voxelSize + 1
    # until that limit all points from the vegetation are allocated to the voxel center

    LAS_ground <<- NA # removing old LAS_veg to ensure performance
    LAS_veg <<- NA # removing old LAS_veg to ensure performance
    gc()

    totalPoints <- allocationCloud@header@PHB$`Number of point records`
    cat(
      thMk(totalPoints), "pts in allocationCloud (blanks) and",
      thMk(beforeLength), "pts in outLAS (voxel seeds)\n"
    )

    dav1 <- Sys.time()
    # two options: assign all Voxel at once,
    #     or if file too big, split it in two and assign voxels in two halves
    #
    if (totalPoints > 500000000) {
      cat("500 Mio points input file is too much!\n -> Reducing the big input file in half:\n")

      maxX <- allocationCloud@header$`Max X`
      minX <- allocationCloud@header$`Min X`
      spanX <- maxX - minX

      maxY <- allocationCloud@header$`Max Y`
      minY <- allocationCloud@header$`Min Y`
      spanY <- maxY - minY

      buffer <- 2 * allocationDistance / 100

      # dividing into lower and upper
      {
        if (spanY > spanX) {
          halfY <- minY + (maxY - minY) / 2

          cat("Lower file: Y from", round(minY, 0), "to incl.", round(halfY, 0), "m ")
          tC1 <- filter_poi(allocationCloud, Y <= halfY)
          as1 <- filter_poi(outLAS, Y <= halfY + buffer)

          upString <- paste("Upper file: Y from excl.", round(halfY, 0), "to", round(maxY, 0), "m ")
          tC2 <- filter_poi(allocationCloud, Y > halfY)
          as2 <- filter_poi(outLAS, Y >= halfY - buffer)
          rm(allocationCloud)
          gc()
        } else {
          halfX <- minX + (maxX - minX) / 2

          cat("Lower file: X from", round(minX, 0), "to incl.", round(halfX, 0), "m ")
          tC1 <- filter_poi(allocationCloud, X <= halfX)
          as1 <- filter_poi(outLAS, X <= halfX + buffer)

          upString <- paste("Upper file: X from excl.", round(halfX, 0), "to", round(maxX, 0), "m ")
          tC2 <- filter_poi(allocationCloud, X > halfX)
          as2 <- filter_poi(outLAS, X >= halfX - buffer)

          rm(allocationCloud)
          gc()
        }
      }


      # Lower
      {
        dt1 <- Sys.time()
        allPoints <- tC1@header@PHB$`Number of point records`
        cat(thMk(allPoints), "pts found.\n")

        {
          cat("Searching neighbors... ")
          te1 <- Sys.time()
          closest <- as.data.table(
            nn2(as1@data[, .(X, Y, Z)],
              query = tC1@data[, .(X, Y, Z)],
              k = 1, searchtype = "radius", radius = (allocationDistance / 100)
            )
          )
          te2 <- Sys.time()
          cat("done in a ")
          print.difftime(round(te2 - te1, 1))
        }


        # add an point identifier for all Million assigned points
        as1@data <- as1@data[, nn.idx := c(1:nrow(as1@data))]

        # assign StemID, runJay and randomCol to every point
        closest <- merge(closest, as1@data[, c("StemID", "runJay", "randomCol", "nn.idx")],
          all.x = TRUE, sort = FALSE
        )

        # remove all points that are too far away
        closest[nn.dists > (allocationDistance / 100), c("StemID", "runJay", "randomCol")] <- 0

        # keep_close_select <- closest$nn.dists <= (allocationDistance/100)
        deletedPoints <- closest[StemID == 0, .N]
        assignedPoints <- (allPoints - deletedPoints)
        nowShare <- round(assignedPoints / allPoints * 100, 2)

        cat(
          "Assigning", thMk(assignedPoints), "pts of total",
          thMk(allPoints), "(nowShare =", nowShare, "%)\n"
        )

        tC1 <- add_lasattribute_manual(tC1, closest$StemID,
          "StemID", "Single Stem ID",
          type = "short"
        )
        tC1 <- add_lasattribute_manual(tC1, closest$runJay,
          "runJay", "number of run added to total",
          type = "short"
        )
        tC1 <- add_lasattribute_manual(tC1, closest$randomCol,
          "randomCol", "Random Stem ID",
          type = "short"
        )

        dt2 <- Sys.time()
        cat("Lower half assigned in a ")
        print.difftime(round(dt2 - dt1, 1))
        cat("\n")
      }

      # Upper
      {
        cat(upString)
        dt1 <- Sys.time()
        allPoints <- tC2@header@PHB$`Number of point records`
        cat(thMk(allPoints), "pts found.\n")

        {
          cat("Searching neighbors... ")
          te1 <- Sys.time()
          closest <- as.data.table(
            nn2(as2@data[, .(X, Y, Z)],
              query = tC2@data[, .(X, Y, Z)],
              k = 1, searchtype = "radius", radius = (allocationDistance / 100)
            )
          )
          te2 <- Sys.time()
          cat("done in a ")
          print.difftime(round(te2 - te1, 1))
        }

        # add an point identifier for all Million assigned points
        as2@data <- as2@data[, nn.idx := c(1:nrow(as2@data))]

        # assign StemID, runJay and randomCol to every point
        closest <- merge(closest, as2@data[, c("StemID", "runJay", "randomCol", "nn.idx")],
          all.x = TRUE, sort = FALSE
        )

        # remove all points that are too far away
        closest[nn.dists > (allocationDistance / 100), c("StemID", "runJay", "randomCol")] <- 0

        # keep_close_select <- closest$nn.dists <= (allocationDistance/100)
        deletedPoints <- closest[StemID == 0, .N]
        assignedPoints <- (allPoints - deletedPoints)
        nowShare <- round(assignedPoints / allPoints * 100, 2)

        cat(
          "Assigning", thMk(assignedPoints), "pts of total",
          thMk(allPoints), "(nowShare =", nowShare, "%)\n"
        )

        tC2 <- add_lasattribute_manual(tC2, closest$StemID,
          "StemID", "Single Stem ID",
          type = "short"
        )
        tC2 <- add_lasattribute_manual(tC2, closest$runJay,
          "runJay", "number of run added to total",
          type = "short"
        )
        tC2 <- add_lasattribute_manual(tC2, closest$randomCol,
          "randomCol", "Random Stem ID",
          type = "short"
        )

        dt2 <- Sys.time()
        cat("Upper half assigned in a ")
        print.difftime(round(dt2 - dt1, 1))
        cat("\n")
      }

      rm(closest)
      gc()

      cat("Combining both halves... ")
      tC1 <- rbind(tC1, tC2)
      allocationCloud <- tC1
      totalPoints <- allocationCloud@header@PHB$`Number of point records`
      cat(thMk(totalPoints), "pts ok.\n")
      rm(tC2, tC1)
      gc()
    } else
    {
      {
        cat("Searching neighbors... ")
        te1 <- Sys.time()
        closest <- as.data.table(
          nn2(outLAS@data[, .(X, Y, Z)],
            query = allocationCloud@data[, .(X, Y, Z)],
            k = 1, searchtype = "radius", radius = (allocationDistance / 100)
          )
        )
        te2 <- Sys.time()
        cat("done in a ")
        print.difftime(round(te2 - te1, 1))
      }

      # add an point identifier for all Million assigned points
      outLAS@data <- outLAS@data[, nn.idx := c(1:nrow(outLAS@data))]

      # assign StemID, runJay and randomCol to every point
      closest <- merge(closest, outLAS@data[, c("StemID", "runJay", "randomCol", "nn.idx")],
        all.x = TRUE, sort = FALSE
      )

      # remove all points that are too far away
      closest[nn.dists > (allocationDistance / 100), c("StemID", "runJay", "randomCol")] <- 0

      # keep_close_select <- closest$nn.dists <= (allocationDistance/100)
      deletedPoints <- closest[StemID == 0, .N]
      assignedPoints <- (totalPoints - deletedPoints)
      nowShare <- round(assignedPoints / totalPoints * 100, 2)

      cat(
        "Assigning", thMk(assignedPoints),
        "pts (nowShare =", nowShare, "%)\n"
      )

      allocationCloud <- add_lasattribute_manual(allocationCloud,
        closest$StemID, "StemID", "Single Stem ID",
        type = "short"
      )
      allocationCloud <- add_lasattribute_manual(allocationCloud,
        closest$runJay, "runJay", "number of run added to total",
        type = "short"
      )
      allocationCloud <- add_lasattribute_manual(allocationCloud,
        closest$randomCol, "randomCol", "Random Stem ID",
        type = "short"
      )
    }
    dav2 <- Sys.time()
    cat("Voxel points assigned in a ")
    print.difftime(round(dav2 - dav1, 1))
    cat("\n")


    # informations, about which seed tree entities have been assigned this round
    nowStems <- unique(allocationCloud$StemID)
    nowStems <- nowStems[nowStems != 0] # 0 refers to the blanks, is no tree at all
    nowStems <- nowStems[order(nowStems)]

    cat("The", length(nowStems), "assigned tree entities ranging from", range(nowStems), "are:\n")
    cat(nowStems)


    blankLAS <- filter_poi(allocationCloud, StemID == 0)
    blankSetLength <- length(blankLAS@data$X)
    allocationCloud <- filter_poi(allocationCloud, StemID != 0)
    cat("\n")


    cat("Writing output files:\n - crowns.laz... ")
    writeLAS(allocationCloud, paste0(crownPath, fileFinder, "_crowns", locationStr, ".laz"))
    cat(thMk(allocationCloud@header@PHB$`Number of point records`), "pts ok.\n")
    outLAS <<- allocationCloud
  } else {
    outLAS_name <<- paste0(crownPath, fileFinder, "_crowns", locationStr, ".laz")
    cat("Creating output files:\n - crowns.laz... ")
    writeLAS(outLAS, paste0(crownPath, fileFinder, "_crowns", locationStr, ".laz"))
    cat(thMk(outLAS@header@PHB$`Number of point records`), "pts ok.\n")
  }

  if (blankSetLength != 0) {
    cat(" - blanks.laz... ")
    writeLAS(blankLAS, paste0(crownPath, fileFinder, "_blanks", locationStr, ".laz"))
    cat(thMk(blankLAS@header@PHB$`Number of point records`), "pts ok.\n")
  }



  gstop <- Sys.time()
  cat("Global ")
  print.difftime(round(gstop - gstart, 1))
  globalTimeDiff <<- paste0(round(gstop - gstart, 1), " ", units(gstop - gstart))
  cat("\n\n\n")
  sink()
  gc()
}
