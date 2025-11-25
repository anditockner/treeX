#' 
#' Runs the function grabDBH( ) for all latest files in one directory
#' 
#' 
#' @param dir_completedInputLists path to the list of exported trees (new trees > 9000)
#' @param dirPath path to all processed files after function extractVegetation( )
#' @param tileClipping global setting for all tiles (standard setting 2x2 for small points, if large forest scan, set to 4x4 or 5x5 tiles)
#' @export
remeasureNewTrees <- function(dir_completedInputLists, appendix = "", fileFinders_selected = "",
                              dirPath = paste0(getwd(), "/"), tileClipping = 2){
  
  {
    library(dplyr)
    
    #plot(1, 1, type = "n", xaxt = "n", yaxt = "n", 
    #     ylab = "", xlab = "", main = "initializing plot area", bty = "n")
    
    #dir_completedInputLists <- "D:/completed_220422_in_all/"
    lists <- list.files(dir_completedInputLists, pattern = "_trees_", full.names = T)
    cat("BATCH REMEASURING NEW TREES\n")
    cat("Analyzing directory with completed lists:\n     ", dir_completedInputLists)
    
    fileFinders_list <- sub("_.*", "", basename(lists))
    
    date_time <- sub(".*_.*_.*_([0-9]{8}_[0-9]{6}).*", "\\1", 
                     basename(lists))
    
    inputFiles <- data.frame(
      fileFinders_list = fileFinders_list,
      date_time = date_time,
      listPath = lists,
      stringsAsFactors = FALSE
    )
    
    inputFiles$date_time <- as.POSIXct(inputFiles$date_time, 
                                       format = "%Y%m%d_%H%M%S")
    
    uniqueFiles_total <- inputFiles %>%
      group_by(fileFinders_list) %>%
      slice_max(order_by = date_time, n = 1) %>%
      ungroup()
    
    head(uniqueFiles_total)
    cat("\ncontaining", nrow(uniqueFiles_total), "unique fileFinders")
    cat(" (of total", length(lists), "lists)\n")
    
    
    if(fileFinders_selected == ""){
      fileFinders_existing <- basename(list.dirs(dirPath, recursive = F))
      fileFinders_existing <- fileFinders_existing[grep("_ALLGO", 
                                                        basename(fileFinders_existing))]
      
      if(length(fileFinders_existing) == 0){
          fileFinders_existing <- basename(list.dirs(dirPath, recursive = F))
          fileFinders_existing <- fileFinders_existing[ !endsWith(fileFinders_existing, "_images") ]
          fileFinders_existing <- fileFinders_existing[ !endsWith(fileFinders_existing, "_in_laz") ]
          fileFinders_existing <- fileFinders_existing[ !endsWith(fileFinders_existing, "_total_ground_veg") ]
          fileFinders_existing <- fileFinders_existing[ !endsWith(fileFinders_existing, "app") ]
          fileFinders_existing <- fileFinders_existing[ !endsWith(fileFinders_existing, "appOld") ]
          fileFinders_existing <- fileFinders_existing[ !endsWith(fileFinders_existing, "parallel_console") ]
          fileFinders_existing
      }
      cat("Working directory", dirPath, "\ncontaining", length(fileFinders_existing), "fileFinders.\n")
      
    } else {
      fileFinders_existing <- fileFinders_selected
      cat("Working in", dirPath, "\nselecting", length(fileFinders_existing), "fileFinders.\n")
    }
    
    
    if(appendix != ""){
      uniqueFiles_total$fileFinders <- paste0(uniqueFiles_total$fileFinders, "_", appendix)
    }
    uniqueFiles <- uniqueFiles_total[
      is.element(toupper(uniqueFiles_total$fileFinders), 
                 toupper(basename(existingFileFinders))), ]
    
    
    if(nrow(uniqueFiles) != length(fileFinders_existing)){
      cat("WARNING - not all lists have a set in the directory!\n")
    }
    
    cat("\n\nProcessing now", nrow(uniqueFiles), "fileFinders:\n")
    cat("Ranging from", uniqueFiles$fileFinders_list[1], "to", 
        uniqueFiles$fileFinders_list[length(uniqueFiles$fileFinders_list)], "\n")
    
  }
  
  for(i in 1:nrow(uniqueFiles)){
    nowFileFinder <- uniqueFiles$fileFinders_list[i]
    nowTreeList <- uniqueFiles$listPath[i]
    
    {
      #LASfile <- NA
      cat("\n###########################\n",
          "#\n",
          "# REMEASURING NEW TREES ", i, "\n", sep = "")
      cat("# FILEFINDER =",nowFileFinder,"\n")
      cat("# COMPLETED LIST =",nowTreeList,"\n")
      cat("# TODAY IS", paste(Sys.time()),"\n")
      cat("#\n")
      cat("#############################\n")
    }
    
    
    try(grabDBH(nowFileFinder, treeList.path = nowTreeList, 
                remeasure = T, tileClipping = tileClipping))
    
  }
  
}




# new idea to prevent merging of close together trees: 
# _______________________________________________________
#
# it is very important that a clean seed set is produced, 
# therefore it would be good to calculate nearest distance to the neighbor, 
# and if that is closer than the enlarged distance, lets say 1.5x radius, 
# then use only radius + a little bit of tolerance, or, better, move away 
# from other tree in opposite direction (to yield a lune halfmoon)
#
#


#' Remeasure best fitting DBH for new trees that have xy and approximate DBH
#'
#'
#' @param treeList.path input path of trees that should be re-measured
#' @param tileClipping how many tiles should be created (default 3 means: 3x3 = 9 tiles), set to 4 (= 16 tiles) for big files
#' @param remeasure creates new cylinderLAS for measuring diameters, set FALSE only if cylinderLAS is already up to date, otherwise some trees might be missing
#' @param allTrees set this TRUE to also measure trees with numbers less than 9000
#' @param new.numbers assign new sequential numbers for circles (1 being the thickest tree)
#' @param frame.up delta z up from 1.30 m for circle fitting
#' @param frame.down delta z down from 1.30 m for circle fitting
#' @param frame.rad re-locating circles: factor of which input DBH will be multiplied to extend circle-fitting
#' @param frame.minD minimum distance (in meter) to fetch the circle around for re-locating, set to 1 m for a tolerant fit (mis-matching systems)
#' @export
grabDBH <- function(fileFinder, treeList.path = NA, 
                    outFolder = getwd(),
                    tileClipping = 3, 
                    overWriteDBHlist = TRUE, # shall old trees_dbh.txt be replaced for segmenting new volumes
                    keepBorderTrees = F, # set TRUE if there are many trees outside of DTM (check if thats a good idea)
                    remeasure = T, 
                    allTrees = F, new.numbers = F, 
                    frame.up = 0.3, frame.down = 0.3, 
                    frame.rad = 1.3, frame.minD = 0.05){
  library(treeX)
  library(mgcv)
  library(spatstat)
  
  
  
  
  if(2==1){
    
    dirPath <- "D:/Species/awayyyyy/"
    treeList.path <- "D:/trees.txt"
    fileFinder <- "KBKI2"
    remeasure = T
    allTrees = F
    new.numbers = F
    tileClipping = 3
    
    frame.up = 0.3
    frame.down = 0.3
    frame.rad = 1.6
    # OLD #frame.rad.impr = 3
    frame.minD = 1
    
    
    
    allTrees = F
    new.numbers = F 
  }
  
  groundPath <- "_total_ground_veg/"
  dbhPath <- paste0(dirPath, fileFinder, "_ALLGO_100to300_dbh/")
  if(!dir.exists(dbhPath)){
    dbhPath <- paste0(dirPath, fileFinder, "/")
  }
  if(!dir.exists(dbhPath)){
    dir.create(dbhPath)
  }
  circlePath <- paste0(dbhPath, "circleFits/")
  if(!dir.exists(circlePath)) dir.create(circlePath)
  
  sink(paste0(dbhPath,"grabDBH_",format(Sys.time(), "%Y%m%d_%H%M"),"_Rcons.txt"), append = TRUE, split = TRUE)
  
  cat("Starting with dbh re-measurement for set",fileFinder,"\n")
  
  
  cat("Today is", format(Sys.time()), "\n")
  
  allStart <- Sys.time()
  # Seed DBH Measurement mode, here the DBH will be calculated
  # from a 80 cm cylinder with frame up and down borders (30 cm pre-setting)
  
  #plot(clustCloud)
  
  #plot(clustCloud@data$X, clustCloud@data$Y, pch = ".", col = clustCloud@data$Intensity, asp = 1)
  
  
  
  if(is.na(treeList.path)){
    treeList.path <- paste0(dbhPath, "trees_dbh.txt")
  }
  
  tryCatch({
    metaList <- read.table(treeList.path, header = T, sep = "\t")
  }, error=function(cond) {
    metaList <<- read.csv(treeList.path)
  })
  
  if(is.element("Xlm", colnames(metaList))){
    # parsing Maissau lists into our format
    metaList$x <- metaList$Xlm
    metaList$y <- metaList$Ylm
    metaList$z <- NA
    metaList$id <- metaList$BaumNummer
    metaList$species <- metaList$FK_Baumart
    metaList$dbh <- metaList$DBH
    # delete columns here
    metaList <- metaList[, !is.element(colnames(metaList), 
                                       c("Xlm", "Ylm", "Zlm", "BaumNummer", "FK_Baumart", "DBH"))]
    
    metaList$speciesNum <- metaList$species
    
  }
  
  cat(length(metaList$id), "trees found in trees_dbh.txt.\n")
  
  if(!allTrees){
    cat("Only remeasuring", sum(metaList$id >= 9000), "trees with an id > 9000.\n")
  }
  
  
  metaList$oldZ <- metaList$z
  
  cat("Fixing the z-coordintes of the newly found trees.\n")
  tryCatch(
    {
      # read in raster file
      dtm_z <- raster(paste0(groundPath, fileFinder, "_ground_min.grd"))
      metaList$z <- round(extract(dtm_z, SpatialPoints(data.frame(x = metaList$x, y = metaList$y))) + 1.3, 3)
      rems <- sum(is.na(metaList$z))
      if(rems != 0){
        plot(metaList$y ~ metaList$x, asp = 1, cex = 0.4, 
             pch = 16, col = "grey", main = paste0(fileFinder, " - rems: ", rems, " trees"))
        points(metaList$y[is.na(metaList$z)] ~ metaList$x[is.na(metaList$z)], 
               cex = 0.4, pch = 16, col = "red")
        if(!keepBorderTrees){
          cat("Removing", rems, "trees that were out of the DTM area.\n")
          metaList <- metaList[!is.na(metaList$z),]
        } else {
          cat("Removing prohibited - ", rems, "border trees set to z = 0.\n")
          metaList$z[is.na(metaList$z)] <- 0
        }
      }
    }, error = function(error_condition) {
      cat("Error in reading the dtm-model!")
      return()
    })
  
  metaList <- metaList[order(metaList$dbh, decreasing = T),]
  # best way to remove columns: 
  #metaList <- metaList[, !is.element(colnames(metaList), 
  #                                   c("height", "vol", "fz"))]
  if(!is.element("species", colnames(metaList))){
    metaList$species <- "NE"
  }
  if(is.element("speciesNum", colnames(metaList))){
    treeSpecies <- function(number) {
      species <- "YY"
      
      species <- switch (paste0(number),
                         "10" = "Bu",
                         "1" = "Fi",
                         "2" = "Ta",
                         "3" = "La",
                         "11" = "Ei",
                         "4" = "Ki",
                         "5" = "Sk",
                         "25" = "Pa",
                         "27" = "We",
                         "30" = "XL",
                         "19" = "Bi",
                         "12" = "Hb",
                         "18" = "Ki",
                         "14" = "Ah",
                         "15" = "Ul",
                         "13" = "Es",
                         "XX")
      
      return(species)
    }
    
    treeSpecies <- Vectorize(treeSpecies)
    
    metaList$species <- treeSpecies(metaList$speciesNum)
    
  }
  
  metaList$species[is.na(metaList$species)] <- "NE"
  metaList$species[metaList$species == "XX"] <- "NE"
  
  table(metaList$species)
  
  if(!is.element("comment", colnames(metaList))){
    metaList$comment <- ""
  }
  table(metaList$comment)
  
  
  # FIXING comment MISSINGNOS
  metaList$comment[is.na(metaList$comment)] <- 0
  if(!is.element("comment", colnames(metaList)) || sum(metaList$comment == 0) > 0){
    cat("We need new random numbers, as there were none in the stem list!\n")
    rnlist <- sample(1:65535, length(metaList$x))
    metaList$comment <- rnlist
    write.table(metaList, file = paste0(dbhPath, fileFinder, "_inTrees.txt"), 
                row.names = FALSE, sep = "\t")
  }
  
  # CREATE NEW SEED LAS SET
  # frame.up <- 0.3
  # frame.down <- 0.3
  # frame.rad <- 1.6
  # frame.rad.impr <- 3 # if tree is smaller than 5 cm, radius factor is x2
  
  
  
  if(!remeasure && file.exists(paste0(dbhPath, "seedLAS_cylinders.las"))){
    cat("Reading old seedLAS_cylinders.las file... ")
    seedLAS <- readLAS(paste0(dbhPath, "seedLAS_cylinders.las"), select = "xyzcit0")
    cat("done!\n")
    
  } 
  else 
  {
    cat("\nCutting a cylinder of each seed tree:\n")
    
    cat("Settings for seedLAS_cylinders: fr.up =", frame.up, "  fr.down =", frame.down, 
        "  fr.rad =", frame.rad," fr.minDist = ", frame.minD, "m.\n")
    
    totalCloud.name <- paste0(dirPath, groundPath, fileFinder, "_raw_veg.laz")
    groundCloud.name <- paste0(dirPath, groundPath, fileFinder, "_ground.laz")
    if(exists("LAS_veg_name") && LAS_veg_name == fileFinder){
      if(exists("LAS_veg") && !is.na(LAS_veg)){
        cat("Using old vegetationCloud!\n")
        totalCloud <- LAS_veg
      } else {
        cat("Reading in totalCloud from", totalCloud.name, "\n")
        totalCloud <- readLAS(totalCloud.name, select = "xyzcit0")
      }
      if (exists("LAS_ground") && !is.na(LAS_ground)) {
        cat("Also using old groundCloud!\n")
        totalCloud <- rbind(totalCloud, LAS_ground)
      } else {
        cat("Reading in groundCloud from", groundCloud.name, 
            "\n")
        groundCloud <- readLAS(groundCloud.name, select = "xyzcit0")
      }
    } else {
      cat("Reading in totalCloud from", totalCloud.name, "\n")
      totalCloud <- readLAS(totalCloud.name, select = "xyzcit0")
      cat("Reading in groundCloud from", groundCloud.name, "\n")
      groundCloud <- readLAS(groundCloud.name, select = "xyzcit0")
      
      totalCloud <- rbind(totalCloud, groundCloud)
      
    }
    
    totalCloud <- add_lasattribute(totalCloud, 0, "StemID", "Single Stem ID")
    totalCloud <- add_lasattribute(totalCloud, 0, "comment", "Random Stem ID")
    
    seedLAS <- filter_poi(totalCloud, StemID != 0) # getting every detected stempoint as seed
    
    
    
    # IDEA: plotting the overall outline of the LAS to remove non-tree points
    
    #plot(metaList$y ~ metaList$x, asp = 1, cex = 0.3)
    #abline(v = c(totalCloud@header@PHB$`Max X`, totalCloud@header@PHB$`Min X`), col = "red")
    #abline(h = c(totalCloud@header@PHB$`Max Y`, totalCloud@header@PHB$`Min Y`), col = "red")
    
    if (!keepBorderTrees) {
      metaList <- metaList[metaList$x >= totalCloud@header@PHB$`Min X` & metaList$x <= totalCloud@header@PHB$`Max X`,]
      metaList <- metaList[metaList$y >= totalCloud@header@PHB$`Min Y` & metaList$y <= totalCloud@header@PHB$`Max Y`,]
      
      cat("Reduced the metaList to", length(metaList$id), "trees in",fileFinder,"\n\n")
      
    }
    
    #for every tree! put into seed set
    i <- 1
    
    
    if(tileClipping == 0){
      cat("Efficiency improvement: Cutting slice from", min(metaList$z, na.rm = TRUE), "to",
          max(metaList$z, na.rm = TRUE), "m.\n")
      cat("Reducing total", thMk(totalCloud@header@PHB$`Number of point records`), "points to ")
      reducedTotalCloud <- filter_poi(totalCloud, Z > min(metaList$z, na.rm = TRUE) - frame.down,
                                      Z < max(metaList$z, na.rm = TRUE) + frame.up)
      cat("cut", thMk(reducedTotalCloud@header@PHB$`Number of point records`), "points.\n")
    }
    
    start <- Sys.time()
    dropoutList <- metaList[0, ]
    metaListSafe <- metaList
    #metaList <- metaListSafe
    
    
    ## TILE GENERATION
    
    tileNumber <- tileClipping #5x5 = 25
    tileBuffer <- 1
    
    xMax <- max(metaList$x) + 0.1
    xMin <- min(metaList$x)-0.1
    yMax <- max(metaList$y) + 0.1
    yMin <- min(metaList$y)-0.1
    
    totalSize.x <- xMax - xMin
    totalSize.y <- yMax - yMin
    
    gridX <- totalSize.x / tileNumber
    gridY <- totalSize.y / tileNumber
    
    
    
    cat("Spanning net of", tileNumber*tileNumber, "tiles x", totalSize.x, 
        "and y", totalSize.y, "m (buffer is", tileBuffer, "m).\n")
    cat(" -> Single tile has x", gridX, "and y", gridY, "m.\n")
    
    xLs <- seq(from = xMin, to = xMax, by = gridX)
    yLs <- seq(from = yMin, to = yMax, by = gridY)
    
    alltexts <- merge(xLs, yLs)
    gridTileName <- function(x, y, noLetters = FALSE){
      ydif <- yMax - y
      xdif <- x - xMin
      
      number <- trunc(xdif / gridX)
      row <- trunc(ydif / gridY)
      outName <- (number + tileNumber*row)
      return(outName)
    }
    
    
    alltexts$tileName <- (gridTileName(alltexts$x + 0.01, alltexts$y-0.01, noLetters = TRUE))
    
    
    #clip(x1 = XU, x2 = XO, y1 = YU, y2 = YO)
    #abline(h = c(YO, YU), v = c(XO, XU), col = "blue")
    
    
    pdf(file = paste0(dbhPath, "tiles", tileClipping, "x", tileClipping, ".pdf"), width = 7, height = 7)
    tryCatch(
      {
        cexFactor <- 5 / mean(gridX + gridY)
        
        plot(metaList$x, metaList$y, type = "n", xlim = c(xMin, xMax), ylim = c(yMin, yMax), asp = 1)
        clip(x1 = xMin, x2 = xMax, y1 = yMin, y2 = yMax)
        abline(h = yLs, v = xLs, col = "grey")
        text(x = alltexts$x + gridX / 2,
             y = alltexts$y-gridY / 2,
             labels = alltexts$tileName, cex = 1.2)
        text(x = metaList$x + 0.5 + cexFactor,
             y = metaList$y-0.5-cexFactor,
             labels = metaList$id, cex = cexFactor*2)
        points(metaList$x, metaList$y, cex = cexFactor)
        
        dev.off()
      }, error=function(e){
        dev.off()
        cat("Error when writing the .pdf, closing device.\n")
      })
    
    
    metaList$tileName <- strtoi(gridTileName(metaList$x, metaList$y, noLetters = TRUE), base = 10)
    cat("Trees per tile: ")
    print(table(metaList$tileName))
    cat("\n")
    tiles <- unique(metaList$tileName)
    tiles <- tiles[order(tiles)]
    
    cat("\nProcessing seeds:\n")
    for(a in 1:length(tiles)){
      if(a %% tileClipping == 0){
        rm(tilesLAS)
        gc()
      }
      cat("-- tile", tiles[a])
      #sta <- Sys.time()
      
      metaList.sub <- metaList[metaList$tileName == tiles[a], ]
      
      xu <- alltexts$x[alltexts$tileName == tiles[a]]
      xo <- xu + gridX
      yo <- alltexts$y[alltexts$tileName == tiles[a]]
      yu <- yo - gridY
      
      tilesLAS <- filter_poi(totalCloud, X > xu - tileBuffer, X < xo + tileBuffer,
                             Y > yu - tileBuffer, Y < yo + tileBuffer)
      cat("", tilesLAS@header@PHB$`Number of point records`, "pts")
      
      
      for(b in 1:length(metaList.sub$id)){
        
        i <- b
        if(b == 1){
          cat(" - ")
        } else {
          cat(", ")
        }
        cat(metaList.sub$id[b])
        #sta <- Sys.time()
        nowClust <- metaList.sub[b, ]
        # print(nowClust)
        nowZ <- nowClust$z
        # cat("b")
        
        
        
        cutRadius <- nowClust$dbh / 200 * frame.rad
        if(cutRadius < frame.minD) { 
          cutRadius <- frame.minD
        }
        
        clustCloud <- clip_circle(tilesLAS, xcenter = nowClust$x, ycenter = nowClust$y,
                                  radius = cutRadius)
        # cat("c", clustCloud@header@PHB$`Number of point records`)
        #cat("Reduced seed to", thMk(clustCloud@header@PHB$`Number of point records`), "points.\n")
        clustCloud <- filter_poi(clustCloud, Z > nowClust$z - frame.down, Z < nowClust$z + frame.up)
        
        if(clustCloud@header@PHB$`Number of point records` == 0){
          cat(" DISCARD: no seed points at DBH found (!)")
          dropoutList[length(dropoutList[, 1]) + 1, ] <- nowClust
          next
        }
        
        # cat("d", clustCloud@header@PHB$`Number of point records`)
        clustCloud@data$StemID <- nowClust$id
        
        # cat("e", clustCloud@header@PHB$`Number of point records`)
        try({clustCloud@data$comment <- nowClust$comment})
        
        # cat(clustCloud@header@PHB$`Number of point records`, "points. - -\n")
        seedLAS@data <- rbind(clustCloud@data, seedLAS@data)
        #sto <- Sys.time()
        #print.difftime(sto-sta)
      }
      cat("\n\n")
    }
    
    write.table(dropoutList, file = paste0(dbhPath, "trees_dropped.txt"),
                row.names = FALSE, sep = "\t")
    
    
    if(keepBorderTrees){
      table(metaList$comment)
      metaList$comment <- ""
      metaList$comment[is.element(metaList$id, dropoutList$id)] <- "skip"
      
      table(metaList$comment)
      cat(nrow(dropoutList), "trees have no points within", frame.minD, "m (comment = \"skip\").\n")
    } else {
      if(length(dropoutList[, 1])!=0){
        cat("\nThese seeds were deleted as false positives, no points present at DBH region:\n")
        print(dropoutList)
      } else {
        cat("\n")
      }
      cleanList <- metaList[!is.element(metaList$id, dropoutList$id),]
      write.table(cleanList, file = paste0(dbhPath, "trees_nonDropped.txt"),
                  row.names = FALSE, sep = "\t")
    }
    
    cat("Generating seedLAS_cylinders.las... \n")
    seedLAS <- add_lasattribute(seedLAS, 1L, "comment", "random Col")
    writeLAS(seedLAS, paste0(dbhPath, "seedLAS_cylinders.las"))
    warning("TODO 04-28-2022: generate seedLAS with new and updated numbers for further segmentation!!!\n")
  }
  
  
  
  #plot(seedLAS)
  idList <- sort(unique(seedLAS@data$StemID))
  clustList <- data.frame("id" = idList, "dbh_ref" = NA)
  
  
  metaList <- merge.data.frame(clustList, metaList, by = "id")
  
  if(allTrees){
    metaList <- metaList[order(metaList$dbh, decreasing = T),]
  }
  
  clustList <- data.frame("id" = metaList$id, "dbh_ref" = metaList$dbh)
  clustList$dbh <- -1
  clustList$dbh.gam <- -1
  clustList$sd.gam <- -1 
  clustList$var.dbh <- -1 
  clustList$var.dbh.gam <- -1 
  clustList$var.sd.gam <- -1      
  clustList$x <- -1
  clustList$x.ref <- metaList$x
  clustList$y <- -1
  clustList$y.ref <- metaList$y
  clustList$slices <- -1
  {
    cat("Distribution of dbhs:")
    clustList$dbhGroup <- ceiling(clustList$dbh_ref/10)*10
    print(table(clustList$dbhGroup))
  }
  clustList$move <- -1
  
  
  if(!allTrees & max(clustList$id) < 9000){
    warning("This set has no trees to be measured. Specify them by id 9000+!")
    if(overWriteDBHlist){
      
      metaList.name <- paste0(dbhPath, "trees_dbh.txt")
      write.table(metaList, file = metaList.name,
                  row.names = FALSE, sep = "\t")
      warning("Only input file list copied as trees.txt!")
    }
    
    sink()
    
    
    return()
  }
  if(allTrees){
    firstTreeMeasured <- 1
  } else {
    firstTreeMeasured <- which(clustList$id >= 9000)[1]
  }
  
  
  n.remeasure <- (length(clustList$id) - firstTreeMeasured+1)
  #clustList[firstTreeMeasured,]
  cat("Re-measuring", n.remeasure, "newly assigned clusters.\n")
  
  
  threshold_percent <- 50
  cat("Intensity filtering per seed is",threshold_percent,"%.\n\n")
  
  threshold_percent_kde <- 94
  cat("KDE filtering per slice is",threshold_percent_kde,"%.\n\n")
  
  #nowId <- 9208
  #i <- which(clustList$id == nowId)
  
  selectorius <- c(1:4, 9, 11, 13, 15)
  cat("\t")
  cat(paste(colnames(clustList)[selectorius]), sep = "\t")
  cat("\n")
  
  
  # measure this all trees
  for(i in firstTreeMeasured:length(clustList$id)){
    
    dbhStart <- Sys.time()
    
    cat(paste0(i - firstTreeMeasured + 1, "/", n.remeasure, "\t"))
    
    nowId <- clustList$id[i]
    clustCloud <- filter_poi(seedLAS, StemID == nowId)
    #plot(clustCloud)
    
    # INTENSITY FILTERING
    threshold <- quantile(clustCloud@data$Intensity, 1 - threshold_percent / 100) # all above are 5 %
    clustInt <- filter_poi(clustCloud, Intensity > threshold)
    #cat(" ", clustInt@header@PHB$`Number of point records`, "pts ")
    
    # WEIGHED KERNEL DENSITY 2D
    #image(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate,
    # col = viridis::viridis(20), xlab = "x", ylab = "y")
    kde_dat <- data.frame("x" = clustInt@data$X, "y" = clustInt@data$Y)
    if(length(kde_dat$x) < 10){
      cat("\n",clustList$id[i], "- too few points in seed, skip this one\n")
      next()
    }
    kde_sample <- ks::kde(x = kde_dat, eval.points = kde_dat)
    #n_cols <- 20
    #quantiles <- quantile(kde_sample$estimate, probs = seq(0, 1, l = n_cols + 1))
    #col <- viridis::viridis(n_cols)[cut(kde_sample$estimate, breaks = quantiles)]
    #points(kde$x, pch = ".", col = col) # The data is returned in $x
    #plot(kde$x, col = col, pch = 19, xlab = "x", ylab = "y")
    clustInt <- add_lasattribute(clustInt, kde_sample$estimate, name = "kde", desc = "2D kernel density")
    
    # KDE filtering for circles density
    #hist(kde_sample$estimate)
    threshold_kde <- quantile(clustInt@data$kde, 1 - threshold_percent_kde / 100) # all above are 5 %
    clustInt <- filter_poi(clustInt, kde > threshold_kde)
    
    #plot(clustInt@data$X, clustInt@data$Y, pch = ".", col = clustInt@data$kde, asp = 1)
    
    sliceHeight <- 0.1 # m
    slices <- seq(from = clustInt@header@PHB$`Min Z`,
                  to = clustInt@header@PHB$`Max Z`, by = sliceHeight)
    
    
    par(mfrow = c(1, 1))
    {
      png(filename = paste0(circlePath, nowId, "_circles.png"), width = 1400, height = 800, type = "cairo")
      par(mfrow = c(3, 4), oma = c(1,1,4,1))
      circles <- NULL
      
      for(j in 1:length(slices)){
        par.circle.1 <- rep(-1, 7)
        slice <- filter_poi(clustInt, Z >= slices[j], Z < slices[j] + sliceHeight)
        #cat("There were", thMk(slice@header@PHB$`Number of point records`), "points found.\n")
        #plot(slice@data$Y, slice@data$X, pch = ".", col = slice@data$Intensity)
        #hist(slice@data$Intensity)
        #threshold_percent <- 100
        #threshold <- quantile(slice@data$Intensity, 1 - threshold_percent / 100) # all above are 5 %
        #lasInt <- filter_poi(slice, Intensity >= threshold)
        #cat("After Intensity filtering there are", thMk(lasInt@header@PHB$`Number of point records`), "points remaining (equals only",
        # round(lasInt@header@PHB$`Number of point records` / slice@header@PHB$`Number of point records`*100, 1), "% of original data).\n")
        lasInt <- slice # no intensity filtering
        
        if(lasInt@header@PHB$`Number of point records` < 5){
          #cat("SLICE",j,"- TOO LITTLE POINTS IN SLICE!\n")
          next()
        }
        
        class(try({
          plot(lasInt@data$X, lasInt@data$Y, pch = 16, col = "black", asp = 1, cex = 0.5)
          points(clustList$x.ref[i], clustList$y.ref[i], col = "grey", cex = 2, pch = 4)
          
          # c2 <- CircleFitByLandau(cbind(lasInt@data$X, lasInt@data$Y), IterMAX = 1000)
          #
          # draw.circle(c2[1], c2[2], c2[3] / 2, border = "red")
          #
          
          nowErr <- FALSE
          
          tryCatch({
            co <- capture.output(
              find.clust.circ <- circMclust(datax=lasInt@data$X,
                                            datay=lasInt@data$Y,
                                            method="const",
                                            # nx=1, ny=1, nr=1, # nx=10, ny=10, nr=25,
                                            #nx=25, ny=25, nr=5,
                                            nx=10, ny=10, nr=5,
                                            minsr=0.01, maxsr=0.5,
                                            brmaxit = 200,
                                            nc=2, # nc=1,
                                            # minsd=0.1, maxsd=0.5,
                                            bw=0.01 # bw=0.05
              ))},
            error = function(error_condition) {
              #cat("Error in the circle finding in slice",j,"!\n")
              nowErr <<- TRUE
            })
          
          
          if(nowErr){
            #cat("NEXT SLICE\n")
            plot(1, type = "n")
            next()
          }
          
          if(!is.null(find.clust.circ)){
            par.circle.1 <- as.numeric(bestMclust(find.clust.circ, 1))
            dbh.circle <- round(2 * par.circle.1[3] * 100, 3)
            x.circle <- par.circle.1[1]
            y.circle <- par.circle.1[2]
            mtext(paste0("x=", round(x.circle, 1), " \ny=", round(y.circle, 1),
                         " \nd=", round(dbh.circle, 1), " "),
                  side = 1, adj = 1, line = -2, col = "red", cex = 1.5)
            
            mtext(paste0(" move ", round(100*sqrt((x.circle - clustList$x.ref[i])^2 + (y.circle - clustList$y.ref[i])^2), 1)),
                  side = 1, adj = 0, line = -2, col = "grey", cex = 1.5)
            draw.circle(x.circle, y.circle, dbh.circle / 200, border = "red")
            
            points(x.circle, y.circle, col = "red", cex = 2, pch = 4)
            lines(x = c(x.circle, clustList$x.ref[i]),
                  y = c(y.circle, clustList$y.ref[i]), col = "grey")
            
            par.circle.1[6] <- NA
            par.circle.1[7] <- NA
          } else {
            circles <- cbind(circles, par.circle.1)
            #cat("SLICE",j,"- IMPLEMENT SOMETHING, CIRCLE DOESN'T WORK!\n")
            next()
          }
          
          
          # ESTIMATE GAM DBH FROM POINTS
          if(!is.null(find.clust.circ)){
            
            #plot(plot.j$X, plot.j$Y, asp=1)
            auffalten <- lasInt@data
            x_null <- auffalten$X - x.circle
            y_null <- auffalten$Y - y.circle
            
            auffalten$dist <- sqrt(x_null**2 + y_null**2)
            auffalten$winkel <- atan2(y_null, x_null)
            
            
            
            
            # shift starting point
            if(1==1){
              
              # detect quarter with max points and then sixteenth with max points, take average as starting point
              ho <- hist(auffalten$winkel, breaks = c(seq(from=-pi, to = pi + 0.1, by = pi/2)), plot = F)
              quarter.from <- ho$breaks[which(max(ho$counts) == ho$counts)[1]]
              quarter.to <- quarter.from + pi/2
              auffalten_quart <- auffalten[auffalten$winkel >= quarter.from & auffalten$winkel <= quarter.to,]
              ho2 <- hist(auffalten_quart$winkel, breaks = c(seq(from=quarter.from, to = quarter.to + 0.1, by = pi/8)), plot = F)
              piece.from <- ho2$breaks[which(max(ho2$counts) == ho2$counts)[1]]
              piece.to <- piece.from + pi/8
              shift_maxDensity <- mean(piece.from, piece.to)
              
              
              
              auffalten$winkel_old <- auffalten$winkel
              #plot(auffalten$winkel_old, auffalten$dist, xlim=c(-pi, pi), ylim=c(0, 0.5), xlab="Winkel", ylab="Distanz")
              #abline(v = shift_maxDensity, col = "green")
              
              auffalten$winkel <- auffalten$winkel - (shift_maxDensity + pi)
              auffalten$winkel[auffalten$winkel < -pi] <- auffalten$winkel[auffalten$winkel < -pi] + 2*pi
              
            }
            
            #range(auffalten$winkel)
            
            nowErr <- FALSE
            tryCatch({
              gam_winkel <- gam(dist ~ s(winkel, bs="cc"), data = auffalten)}, 
              error = function(error_condition) {
                #cat("Error in the gam in slice",j,"!\n")
                nowErr <<- TRUE
              })
            
            
            if(nowErr){
              #cat("NEEEXTice!")
              plot(1, type = "n")
              circles <- cbind(circles, par.circle.1)
              next()
            } 
            
            #gam.check(gam_winkel)
            #summary(gam_winkel)
            res_gam <- residuals(gam_winkel) #Residuen vom GAM
            #hist(res_gam)
            sd(res_gam) #Standardabweichung Residuen #0.02072122
            
            #st_res_gam <- (res_gam - mean(res_gam))/sd(res_gam) #Standardisierte Residuen vom Gam
            #hist(st_res_gam)
            
            
            gamfunc <- function(winkel) { predict(gam_winkel, newdata=data.frame("winkel"=winkel)) }
            
            werte <- data.frame("winkel"=seq(-pi, pi, 2*pi/360))
            werte$ergeb <- gamfunc(werte$winkel)
            
            #title(main = paste(cluster.vec[i], "_", cluster2[clust2], "_", cluster3[clust3], ", ", u.grenzen.vec[j], sep=""))
            
            #?gam
            
            vorhersage <- data.frame("winkel"=seq(-pi, pi, 2*pi/360))
            vorhersage$dist <- gamfunc(vorhersage$winkel)
            vorhersage$winkel <- vorhersage$winkel + (shift_maxDensity + pi)
            vorhersage$X <- x.circle + cos(vorhersage$winkel) * vorhersage$dist
            vorhersage$Y <- y.circle + sin(vorhersage$winkel) * vorhersage$dist
            
            #plot(auffalten$X, auffalten$Y, asp=1) #Zeichnen der Punkte
            
            w <- owin(poly=list(x=c(vorhersage$X), y=c(vorhersage$Y)))
            
            #d aus flaeche
            area <- area.owin(w)
            d.gam <- sqrt(area/pi*4) * 100 #d aus flaeche
            
            
            #points(vorhersage$X, vorhersage$Y, col=2, cex=1.5, pch=18)
            # green line is the starting point of the gam
            lines(x = c(x.circle, x.circle + d.gam/200 * cos(shift_maxDensity)),
                  y = c(y.circle, y.circle + d.gam/200 * sin(shift_maxDensity)), col = "green")
            plot(w, add=T, border = "green")
            
            plot(auffalten$winkel, auffalten$dist, xlim=c(-pi, pi), ylim=c(0, 0.5), xlab="Winkel", ylab="Distanz")
            lines(werte$winkel, werte$ergeb, col="green", lwd=2)
            mtext(paste0("sd=", round(sd(res_gam)*100,1),
                         " \nd.gam=", round(d.gam, 1), " "),
                  side = 3, adj = 1, line = -5, col = "green", cex = 1.5)
            
            #title(main = paste(cluster.vec[i], "_", cluster2[clust2], "_", cluster3[clust3], ", ", u.grenzen.vec[j], ", \n d.circ: ",
            #                   round(dbh.circle, 2), ", d.circ2: ", round(dbh.circle.2, 2), ", d.ell: ", round(dbh.ellipse, 2), ", d.gam: ", round(d.gam, 2), sep=""))
            
            d.gam
            par.circle.1[6] <- d.gam
            par.circle.1[7] <- sd(res_gam)*100
            
            
          }
          
          circles <- cbind(circles, par.circle.1)
        }))
        
        
        
        
        
        
        
        
      }
      
      try(mtext(
        paste0(clustList[i,]$id, 
               " dbhs: circ ", round(median(circles[3, ])*200, 1),
               " (vr", round(sd(circles[3, ]*200, na.rm = T),2), ")",
               " gam ", round(median(circles[6, ], na.rm = T), 1), 
               " (vr", round(sd(circles[6, ], na.rm = T),2), ")",
               " sd ", round(median(circles[7, ], na.rm = T), 1), 
               " (vr", round(sd(circles[7, ], na.rm = T),2), ")",
               " --- alt ", round(clustList[i,]$dbh_ref, 1))
        , 
        outer = T, cex = 3
      )
      )
      
      dev.off()
      if(!is.null(circles)){
        x.mean <- round(median(circles[1, ]), 3)
        y.mean <- round(median(circles[2, ]), 3)
        dbh.mean <- round(median(circles[3, ], na.rm = T)*200, 3)
        dbh.gam.mean <- round(median(circles[6, ], na.rm = T), 3)
        gam.sd.mean <- round(median(circles[7, ], na.rm = T), 3)
        dbhStop <- Sys.time()
        timeNB <- as.difftime(dbhStop - dbhStart)
        #cat("DBH circle", dbh.mean, "cm and gam",dbh.gam.mean,"cm in", round(timeNB, 1), units(timeNB), "\n")
        
        clustList[i, ]$dbh <- dbh.mean
        clustList[i, ]$dbh.gam <- dbh.gam.mean
        clustList[i, ]$sd.gam <- gam.sd.mean
        clustList[i, ]$var.dbh <- round(sd(circles[3, ]*200, na.rm = T),2)
        clustList[i, ]$var.dbh.gam <- round(sd(circles[6, ], na.rm = T),2)
        clustList[i, ]$var.sd.gam <- round(sd(circles[7, ], na.rm = T),4)
        # HUOM: Not real variance, but sd just to prevent name confusion with gam sd (scattering)
        clustList[i, ]$x <- x.mean
        clustList[i, ]$y <- y.mean
        clustList[i, ]$slices <- length(slices)
        clustList[i, ]$move <- round(100*sqrt((x.mean - clustList$x.ref[i])^2 + (y.mean - clustList$y.ref[i])^2), 1)
      } else {
        timeNB <- as.difftime(dbhStop - dbhStart)
        #cat("No DBH estimated in", round(timeNB, 1), units(timeNB), "\n")
      }
      
      
      
      
      #dev.off()
      
    }
    cat(paste(clustList[i, selectorius]), sep = "\t")
    cat("\n")
  }
  
  write.table(clustList, paste0(dbhPath, "trees_remeasured_dbh_circleFits.txt"), row.names = F, sep = "\t", dec = ".")
  
  
  
  
  
  # rulebook ####
  {
    cat("\nDeciding now, which new diameters to take as valid.\n")
    # check which new measured to take into
    # CURRENT PROJECT RULE BOOK (ladiwaldi): 
    #    > 140 cm -> GAM dbh is considered invalid 
    #    take GAM 
    #    Label as BAD if: 
    #      REF - GAM > 5 cm absolute
    #         still good if: CIRC / GAM between 0.9 and 1.1 (90 - 110 % relative)
    #         still good if: CIRC - GAM < 5 cm absolute
    #    take CIRC if: 
    #      BAD and REF - CIRC <= 5 cm
    #    DBH < 1 cm = 1 cm
    # # # # # 
    # shift if: 
    #    move distance < radius 
    #    move distance < 20 cm
    
    if(allTrees){
      improveList <- clustList
    } else {
      improveList <- clustList[clustList$id >= 9000,] 
    }
    improveList$dbh.gam[improveList$dbh.gam > 140] <- -11
    improveList$dbh.gam[improveList$dbh.gam < 0] <- -11
    improveList$dbh.gam[is.na(improveList$dbh.gam)] <- -11
    improveList$dbh[improveList$dbh < 0] <- -11
    improveList$dbh[is.na(improveList$dbh)] <- -11
    cat("In total there are", length(improveList$id), 
        "new measured dbhs ( na.gam =",sum(improveList$dbh.gam == -11),
        "and na.circ =", sum(improveList$dbh == -11),").\n")
    
    improveList$diff.dist <- sqrt((improveList$x.ref - improveList$x)^2 + (improveList$y.ref - improveList$y)^2 )
    improveList$diff.dbh <- abs(improveList$dbh_ref - improveList$dbh)
    improveList$diff.dbh.gam <- abs(improveList$dbh_ref - improveList$dbh.gam)
    
    cat("5 cm accepted tolerance to guessed DBH: ")
    acceptToleranceCM <- 5
    #table(improveList$dbh.gam != -11) 
    improveList$dbh.check <- "Fi" #good ones
    improveList$dbh.sugg <- improveList$dbh.gam
    improveList$dbh.bad <- " "
    improveList$dbh.bad[improveList$diff.dbh.gam > acceptToleranceCM] <- "X"
    #improveList$dbh.bad[improveList$dbh.gam == -11] <- "X"
    table(improveList$dbh.bad[improveList$dbh.gam != -11]) 
    
    cat(sum(improveList$dbh.bad[improveList$dbh.gam != -11] != "X"), "assigned\n") 
    table(improveList$dbh.bad) 
    
    
    
    cat("90-110 % accepted ratio to guessed DBH: ")
    #set 90 - 110 % circle/gam equivalence of both measures as acceptable (Buche)
    improveList$dbh.check[improveList$dbh.bad == "X"] <- "Bu" #65 with 90-110% equivalence
    improveList$dbh.est.ratio <- improveList$dbh/improveList$dbh.gam*100
    improveList$dbh.bad[improveList$dbh.est.ratio >= 90 & 
                          improveList$dbh.est.ratio <= 110 & 
                          improveList$dbh.gam != -11] <- " "
    table(improveList$dbh.bad[improveList$dbh.gam != -11])
    cat(sum(improveList$dbh.bad[improveList$dbh.gam != -11] != "X"), "assigned\n") 
    table(improveList$dbh.bad)
    
    
    
    
    cat("5 cm allowed discrepancy to circle DBH: ")
    #set absolute difference 5 cm as acceptable (Ulme)
    improveList$dbh.check[improveList$dbh.bad == "X"] <- "Ul" #6 with 5cm difference
    improveList$dbh.est.ratio.diff <- abs(improveList$dbh - improveList$dbh.gam)
    improveList$dbh.bad[improveList$dbh.est.ratio.diff <= acceptToleranceCM & 
                          improveList$dbh.gam != -11] <- " "
    table(improveList$dbh.bad[improveList$dbh.gam != -11]) 
    #a_Rosalia: 1 bad ones + 133
    #b_Altmuenster: 2 bad ones + 58
    cat(sum(improveList$dbh.bad[improveList$dbh.gam != -11] != "X"), "assigned\n") 
    table(improveList$dbh.bad)
    
    
    
    cat("CIRCLE 5 cm discrepancy to guessed DBH: ")
    # 133 bad gam ones: take the dbh of circular measurements if difference is less than 5 cm
    improveList$dbh.check[improveList$dbh.bad == "X"] <- "Ta" #51 dbh circle values
    choice.bad.takeCircledbh <- improveList$dbh.bad == "X" & improveList$diff.dbh <= 5
    improveList$dbh.bad[choice.bad.takeCircledbh] <- " "
    table(improveList$dbh.bad) #82 left
    improveList$dbh.sugg[choice.bad.takeCircledbh] <- improveList$dbh[choice.bad.takeCircledbh]
    table(improveList$dbh.bad) #82 no good measurement, Altmuenster: 18 no good
    cat(sum(improveList$dbh.bad != "X"), "assigned\n") 
    #improveList[improveList$dbh.bad=="X",]
    
    
    
    
    # trees smaller than 1 cm dbh
    #
    thoseSmallerThanOneCm <- improveList$dbh.sugg != -11 & improveList$dbh.sugg < 1
    if(sum(thoseSmallerThanOneCm) > 0){
      cat("Trees smaller than 1 cm correction for",sum(thoseSmallerThanOneCm),"trees.")
      improveList[improveList$dbh.sugg != -11 & improveList$dbh.sugg < 1,]
      improveList$dbh.sugg[improveList$dbh.sugg != -11 & improveList$dbh.sugg < 1] <- 1 # set to 1 cm
      
    }
    improveList$dbh.check[improveList$dbh.sugg == -11] <- "Ks" #82 not measured, stay estimates
    print(table(improveList$dbh.check))
    sum(improveList$dbh.gam == -11) #133 no gam
    sum(improveList$dbh.sugg == -11) #82 not measured (51 from circle dbh)
    notImprovedList <- improveList[improveList$dbh.sugg == -11,]
    #table(notImprovedList$dbh_ref) # only up to 5 cm dbh, Altmuenster one very outside 30 cm (not anymore in area)
    
    
    
    cat("There are",length(notImprovedList$dbh),"new trees without improved dbh.\n")
    cat(notImprovedList$id, sep = " - ")
    cat("\n")
    
    
    #improveList <- improveList[improveList$dbh.gam != -1,]
    #improveList$id[improveList$dbh.check == "X"] # check these ones, they have
    #improveList[improveList$dbh.sugg > 100,]
    
    
    
    
    
    improveList$shift.x <- NA
    improveList$shift.y <- NA
    
    # SHIFT ALL
    shiftAll <- FALSE
    if(shiftAll){
      cat("Force to shift all trees!\n")
      shift.those <- improveList$move/100 < frame.minD
      improveList[!shift.those,]
      cat("\nShifting",sum(shift.those, na.rm = T),"to new center, at most", frame.minD * 100, "cm.")
      
    } else {
      shift.those <- improveList$diff.dist < improveList$dbh_ref/100*2
      shift.those.also <- improveList$diff.dist < 0.2
      shift.those <- shift.those | shift.those.also
      table(shift.those)
      cat("\nShifting",sum(shift.those, na.rm = T),"to new center, at most 20 cm or 2x dbh.")
      
      
      
      # check the ones that need to be manually shifted
      check.shift.man <- improveList[!is.na(improveList$diff.dist) & improveList$diff.dist > 0.2,]
      check.shift.man <- check.shift.man[order(check.shift.man$id),]
      
      check.shift.man[,c(1:3,8,9:12,14)]
      
      #shift.man <- c(488, 9778, )
    }
    
    
    
    clustList$x[clustList$x == -1] <- NA
    clustList$y[clustList$y == -1] <- NA
    
    unshifted <- !shift.those & improveList$dbh.sugg != -11
    if(sum(unshifted > 0)){
      cat("Not shifted trees of improved file: \n")
      print(improveList[unshifted,c(1,2,3,4,9:12,14)])
    } else {
      cat("All trees shifted!")
    }
    
    improveList$shift.x[shift.those] <- improveList$x[shift.those]
    improveList$shift.y[shift.those] <- improveList$y[shift.those]
    improveList$shift.col <- "Ta"
    improveList$shift.col[shift.those] <- "Vb"
    
    mergeImprove <- data.frame("id" = improveList$id, 
                               "dbh.ref" = improveList$dbh_ref, 
                               "dbh.sugg" = improveList$dbh.sugg, 
                               "dbh.check" = improveList$dbh.check, 
                               "shift.x" = improveList$shift.x, 
                               "shift.y" = improveList$shift.y, 
                               "shift.col" = improveList$shift.col)
    mergeImprove$dbh.sugg[mergeImprove$dbh.sugg == -11] <- NA
    write.table(mergeImprove, paste0(dbhPath, "trees_rem_shift_compare.txt"), row.names = F, sep = "\t", dec = ".")
    
    
    
    
    
    
    
    # applying new measurements to trees_dbh.txt
    
    cat("\nCombining into new trees_dbh.txt...\n")
    trees_out <- merge.data.frame(metaList, mergeImprove, by = "id", all.x = T)
    
    onlyShift <- TRUE
    if(!onlyShift){
      trees_out$dbh.old <- trees_out$dbh
      trees_out$dbh[!is.na(trees_out$dbh.sugg)] <- trees_out$dbh.sugg[!is.na(trees_out$dbh.sugg)]
      table(trees_out$dbh.check)
      trees_out$dbh.check[is.na(trees_out$dbh.check)] <- "Ei" # old kept dbh of CG found trees
      trees_out$shift.col[is.na(trees_out$shift.col)] <- "Ei" # old not shifted of CG found trees
      trees_out <- trees_out[,-2]
    }
    shift.those <- !is.na(trees_out$shift.x)
    table(shift.those)
    
    trees_out$x_old <- trees_out$x
    trees_out$y_old <- trees_out$y
    
    trees_out$x[shift.those] <- trees_out$shift.x[shift.those]
    trees_out$y[shift.those] <- trees_out$shift.y[shift.those]
    
    trees_out$dbh[!is.na(trees_out$dbh.sugg)] <- round(trees_out$dbh.sugg[!is.na(trees_out$dbh.sugg)],1)
    
    
    
    
    
    trees_out$oldID <- trees_out$id
    if(new.numbers){
      
      trees_out <- trees_out[order(trees_out$dbh, decreasing = T),]
      # #trees_out$oldID <- metaList$id
      # trees_out$newID <- c(1:length(trees_out$id))
      trees_out$id <- c(1:length(trees_out$id))
      trees_out <- trees_out[order(trees_out$id, decreasing = F),]
      write.table(trees_out, paste0(dbhPath, "trees_remeasured_total.txt"), row.names = F, sep = "\t", dec = ".")
      
      cat(" -> file", paste0("\"trees_remeasured_total.txt\""), "successfully created!\n")
      
      
      
    }
    
    
    #trees_out$species <- trees_out$dbh.check
    #trees_out$oldID <- trees_out$id
    #trees_out$species <- metaList$species
    #trees_out$id <- metaList$id
    #trees_out <- trees_out[order(trees_out$id, decreasing = F),]
    
    
    appList <- trees_out[,c(which(colnames(trees_out) == "x"),
                            which(colnames(trees_out) == "y"),
                            which(colnames(trees_out) == "z"),
                            which(colnames(trees_out) == "id"),
                            which(colnames(trees_out) == "comment"),
                            which(colnames(trees_out) == "species"),
                            which(colnames(trees_out) == "dbh"),
                            which(colnames(trees_out) == "height"),
                            which(colnames(trees_out) == "vol"),
                            which(colnames(trees_out) == "fz"),
                            which(colnames(trees_out) == "inside"),
                            which(colnames(trees_out) == "selected"))]
    
    
    appList.dropped <- dropoutList[,c(which(colnames(dropoutList) == "x"),
                                      which(colnames(dropoutList) == "y"),
                                      which(colnames(dropoutList) == "z"),
                                      which(colnames(dropoutList) == "id"),
                                      which(colnames(dropoutList) == "comment"),
                                      which(colnames(dropoutList) == "species"),
                                      which(colnames(dropoutList) == "dbh"),
                                      which(colnames(dropoutList) == "height"),
                                      which(colnames(dropoutList) == "vol"),
                                      which(colnames(dropoutList) == "fz"),
                                      which(colnames(dropoutList) == "inside"),
                                      which(colnames(dropoutList) == "selected"))]
    
    appList <- rbind(appList, appList.dropped)
    
    appList <- appList[order(appList$id, decreasing = F),]
    
    write.table(appList, 
                file = paste0(outFolder, "/", 
                              tolower(strRep(fileFinder, "_", "+")), 
                              "_trees_all_", format(Sys.time(), "%Y%m%d_%H%M%S"), 
                              "_remeasured.txt"),
                row.names = FALSE, sep = "\t")
    
    
    
    if(overWriteDBHlist){
      
      
      outList <- trees_out[,c(which(colnames(trees_out) == "x"),
                              which(colnames(trees_out) == "y"),
                              which(colnames(trees_out) == "z"),
                              which(colnames(trees_out) == "id"),
                              which(colnames(trees_out) == "comment"),
                              which(colnames(trees_out) == "dbh"),
                              which(colnames(trees_out) == "dbh.ref"),
                              which(colnames(trees_out) == "oldID"),
                              which(colnames(trees_out) == "species"),
                              which(colnames(trees_out) == "dist"),
                              which(colnames(trees_out) == "angle"),
                              which(colnames(trees_out) == "inside"))
      ]
      outList$x <- round(outList$x, 3)
      outList$y <- round(outList$y, 3)
      outList$z <- round(outList$z, 3)
      outList$dbh <- round(outList$dbh, 1)
      outList$dbh.ref[is.na(outList$dbh.ref)] <- outList$dbh[is.na(outList$dbh.ref)]
      
      metaList.name <- paste0(dbhPath, "trees_dbh.txt")
      metaList.name.old <- paste0(dbhPath, "trees_dbh_uncorrected_", format(Sys.time(), "%Y%m%d_%H%M"), ".txt")
      if(file.rename(from = metaList.name, to = metaList.name.old)){
        cat("Old file renamed to", basename(metaList.name.old), "!\n")
      } else {
        cat("Something went wrong with renaming sophisticated trees_dbh.txt (file was not present)...\n")
      }
      
      cat("File", (metaList.name), "succesfully created.\n")
      write.table(outList, file = metaList.name,
                  row.names = FALSE, sep = "\t")
    }
    
  }
  
  
  
  
  
  cat("\nDBH re-measurement for set",fileFinder,"done in a ")
  print.difftime(round(Sys.time() - allStart, 1))
  sink()
  cat("\n\n\n")
}


