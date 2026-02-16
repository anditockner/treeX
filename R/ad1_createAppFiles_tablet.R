



#' Creation of APP files (background image and tree list)
#'
#' Reads in the las slices of the given set in fileFinder, 
#' and uses trees_dbh.txt or trees_measured.txt to export the tree list
#'
#' @param fileFinder user defined name to find this certain dataset in further processing
#' @param changeFileFinder rename app background to this name
#' @param quickPreview set TRUE to plot only reduced density of slice
#' @param drawGround set FALSE to omitt the voxelized ground points lower than 1 m
#' @param fullGroundGrey set TRUE to plot full ground points from 1 m lower in grey (instead of voxelized points)
#' 
#'
#' @param pixelUnit_cm how many cms are represented in one pixel - default 1 (100 pixel equals one meter), also 2 works for smaller plots
#' @param circleRadius radius of the additional circle that is plot in blue on top of the background file - recommended: 5 m more than desired sample plot radius
#' @param fixedLimit limit in meters that defines boundary of background file - for round circles set it to the clipped file size, so that all backgrounds have the same dimensions
#' @param exportClippedLAS if set TRUE, then a laz file will be exported according to clip.radius
#' @export
createAppFiles <- function(fileFinder = NA, 
                           drawGround = T, 
                           fullGroundGrey = F,
                           drawRedSlice = T, 
                           pixelUnit_cm = 2, 
                           
                           createBGR_pic = TRUE, 
                           createJPG = TRUE, jpgQuality = 0.3, # 0.3 means 30%, for smaller files can also increase to 50%
                           createTIFF = FALSE, 
                           dirPath = paste0(getwd(), "/"),
                           changeFileFinder = "", addTrajLAZ = F, 
                           drawTraj = FALSE, 
                           greySpots = data.frame(),
                           drawLines = data.frame(), 
                           cutWindow = c(-1000, -1000, 2000, 2000), 
                           reduceName = FALSE, # take only first lettter of long_name_now = lnn
                           onlySuffixFileFinder = FALSE, # take only the last underscore example: project_c2b -> c2b for non-registered point-clouds

                           slices = c("_clusterSlice_100to300.laz"),
                           laz.path = NA, slice.cut.low = NA, slice.cut.high = NA,
                           colorVersion = F, 
                           eraseSpecies = F, highlightNew = F, setSpecies = "", 
                           dtm.path = "", trees.path = "", 
                           isoLines = 0, thin = FALSE,
                           modeApp = 1, selector = "xyzi",
                           writeColoredLAZ = FALSE,
                           circleRadius = 0, fixedLimit = 0,
                           wait5s = FALSE, 
                           quickPreview = F, greenTreeLocations = F){
  
  
  if(1==2){
    # debug
    fileFinder = NA
    changeFileFinder = ""
    reduceName = FALSE # take only first lettter of long_name_now = lnn
    onlySuffixFileFinder = FALSE # take only the last underscore example: project_c2b -> c2b for non-registered point-clouds
    
    slices = c("_clusterSlice_100to300.laz")
    laz.path = NA
    slice.cut.low = NA
    slice.cut.high = NA
    drawGround = T
    colorVersion = F
    eraseSpecies = F
    
    dtm.path = ""
    trees.path = ""
    jpgQuality = 0.3 # 0.3 means 30%, for smaller files can also increase to 50%
    drawRedSlice = T
    thin = FALSE
    
    pixelUnit_cm = 1
    modeApp = 1
    writeColoredLAZ = FALSE
    circleRadius = 0
    fixedLimit = 0
    quickPreview = F
    createBGR_pic = T
    greenTreeLocations = F
  }
  
  if(createTIFF) createJPG <- FALSE
  # for borders
  library(sf)
  library(spatstat)
  library(alphahull)
  
  # for converting png to jpg (faster than writing jpg)
  # correct jpeg
  library(png)
  library(jpeg)
  library(pracma) # for strRep()
  
  writeLASuncompressed <- FALSE 
  # set to TRUE if you prefer .las files
  
  if(!is.na(fileFinder[1])){
    if(length(fileFinder) == 1){
      if(createBGR_pic) {
        if(createTIFF){
          cat("Creating app files (.tiff and .txt) for set", fileFinder, "\n")
        } else if(createJPG){
          cat("Creating app files (.jpg and .txt) for set", fileFinder, "\n")
        } else{
          cat("Creating app files (.png and .txt) for set", fileFinder, "\n")
        }
      } else {
        cat("Creating app files (only .txt) for set", fileFinder, "\n")
      }
      
      outName <- tolower(fileFinder)
      if(changeFileFinder != ""){
        outName <- changeFileFinder
        cat("Changing set name to", outName, "\n")
      } else if(onlySuffixFileFinder){
        parts <- strsplit(outName, "_")
        outName <- parts[[1]][length(parts[[1]])]
        cat("Keeping only last fileFinder as set name ", outName, "\n")
      } else if(reduceName){
        # reduce the name to first letter of every _ segmented word
        # just in case I ever need it...
        parts <- strsplit(outName, "_")
        outName <- unlist(lapply(parts, substr, 1, 1))
        outName <- paste0(stringi, collapse = "")
        cat("Reduce set name to ", outName, "\n")
      }
      
      
      
      
      
    } else {
      if(createBGR_pic) {
          if(createTIFF){
            cat("Creating app files (.tiff and .txt) for", length(fileFinder), "sets:\n")
          } else if(createJPG){
            cat("Creating app files (.jpg and .txt) for", length(fileFinder), "sets:\n")
          } else{
            cat("Creating app files (.png and .txt) for", length(fileFinder), "sets:\n")
          }
      } else {
        cat("Creating app files (only .txt) for", length(fileFinder), "sets:\n")
      }
      cat("        ", paste0(fileFinder, collapse = ",  "))
      
      
      outName <- tolower(paste0(fileFinder, collapse = ""))
      if(changeFileFinder != ""){
        outName <- changeFileFinder
      } else if(reduceName){
        # reduce the name to first letter of every _ segmented word
        # just in case I ever need it...
        outName <- tolower(paste0(fileFinder, collapse = "_"))
        parts <- strsplit(outName, "_")
        outName <- unlist(lapply(parts, substr, 1, 1))
        outName <- paste0(outName, collapse = "")
      }
      
      cat("   ->   merge to \"", outName, "\"\n", sep = "")
    }
    cat("Today is", format(Sys.time()), "\n")
    
    
    
    
    
    
    if(!exists("dirPath")){
      cat("ERROR: NO INPUT PATH FOUND!\nPlease specify the directory in global variable dirPath!\n\n")
      return()
    }
    groundPath <- paste0(dirPath, "_total_ground_veg/")
    
    
    if(drawTraj){
      cat("Searching for trajectory... ")
      trajFiles <- paste0(groundPath, fileFinder, "_traj.txt")
      cat(sum(file.exists(trajFiles)), "found!\n")
      trajFiles <- trajFiles[file.exists(trajFiles)]
      
      allTraj <- data.frame()
      for(i in 1:length(trajFiles)){
        if(!file.exists(trajFiles[i])) next()
        try({
          traj <- read.csv(trajFiles[i], sep = " ")
        })
        traj <- traj[seq(from = 1, to = nrow(traj), by = 25), 
        ]
        traj$col <- rainbow(length(traj[, 1]), end = 0.7, rev = T)
        allTraj <- rbind(allTraj, traj)
        
      }
      plot(allTraj$y ~ allTraj$x, col = allTraj$col, cex = 0.6, pch = 16, asp = 1)
      
      rm(traj)
      gc()
    }
    
    
    #slices <- c("_clusterSlice_100to300.laz", "_clusterSlice_100to200.laz",
    #            "_clusterSlice_300to500.laz", "_clusterSlice_50to200.laz") # print all a slice bgr pic
    #slices <- c("_clusterSlice_100to300.laz", "_clusterSlice_50to200.laz")
    # slices <- c("_clusterSlice_100to300.laz")
    slices_gr <- c("_ground_100cm_vox.laz")
    if(fullGroundGrey) slices_gr <- "_ground.laz"
    slices_grTRAJ <- c("_ground_100cm_vxtraj.laz")
    
   
    
  }   else {
    if(is.na(laz.path)){
      warning("No fileFinder and no laz.path specified!\n")
      return()
    }
    
    
    if(!file.exists(dtm.path)){
      if(dtm.path == ""){
        warning(paste0("No DTM specified - please use parameter dtm.path!\n"))
      } else {
        warning(paste0("Specified DTM is missing! No file in ", dtm.path, "!\n"))
      }
      return()
    } else {
      groundModel.path <- dtm.path
      
    }
    
    
    
    
    tempFileFinder <- strsplit(basename(dtm.path), "_")[[1]][1]
    outName <- tolower(tempFileFinder)
    if(changeFileFinder != ""){
      outName <- changeFileFinder
      cat("Changing set name to", outName, "\n")
    } else if(onlySuffixFileFinder){
      parts <- strsplit(outName, "_")
      outName <- parts[[1]][length(parts[[1]])]
      cat("Keeping only last fileFinder as set name ", outName, "\n")
    } else if(reduceName){
      # reduce the name to first letter of every _ segmented word
      # just in case I ever need it...
      parts <- strsplit(outName, "_")
      outName <- unlist(lapply(parts, substr, 1, 1))
      outName <- paste0(stringi, collapse = "")
      cat("Reduce set name to ", outName, "\n")
    }
    
    if(createBGR_pic) {
      if(createTIFF){
        cat("Creating app files (.tiff and .txt) for", length(laz.path), "sets:\n")
      } else if(createJPG){
        cat("Creating app files (.jpg and .txt) for", length(laz.path), "sets:\n")
      } else{
        cat("Creating app files (.png and .txt) for", length(laz.path), "sets:\n")
      }
    } else {
      cat("Creating app files (only .txt) for", length(laz.path), "files: \n  ")
    }
    cat(basename(laz.path), sep = "\n  ")
    cat("Today is", format(Sys.time()), "\n")
    cat("Using DTM from:", dtm.path, "\n")
    dirPath <- paste0(dirname(laz.path[1]), "/")
    groundPath <- dirPath
    workingPath <- dirPath
    cat("Working folder is:", dirPath, "\n")
    cat("\n")
  }
  
  
  ### PLOTTING BGR IMAGE
  if(dtm.path == ""){
    groundModel.path <- paste0(groundPath,fileFinder,"_ground_min.grd")
  } else {
    groundModel.path <- dtm.path
  }
  #fineModel.path <- paste0(groundPath,fileFinder,"_ground_fine.grd")
  #canopyModel.path <- paste0(groundPath,fileFinder,"_dem.grd")
  if(sum(file.exists(groundModel.path))!=length(groundModel.path)){
    warning(paste("DTM file", groundModel.path, "doesn't exist!\n"))
    #return()
  }
  
  
  if(createBGR_pic){
    
    if(colorVersion){
      cat("\nDrawing colorful image")
    } else {
      if(drawRedSlice){
        cat("\nDrawing red DBH slice")
      } else {
        cat("\nDrawing all black")
        
      }
    }
    if(fullGroundGrey){
      drawGround <- T
      cat(", with full grey ground")
    } else {
      if(drawGround){
        cat(", with ground")
      } else {
        cat(" and no ground")
      }
    }
  } else {
    cat("Background .png is skipped")
  }
  if(eraseSpecies){
    cat(", erasing comment (old randomCol) and species")
  }
  if(setSpecies != ""){
    cat(paste0(", setting all to ", "\"", setSpecies, "\""))
  } else if(highlightNew){
    cat(", >9000 trees are assigned to GX")
  }
  if(writeColoredLAZ){
    cat(", incl. colored ")
    if(writeLASuncompressed){
      cat("LAS")
    } else {
      cat("LAZ")
    }
    cat("-files")
  }
  cat(".\n")
  if (sum(cutWindow == c(-1000, -1000, 2000, 2000)) != 4) {
    cat("Cutting to a window from x =", cutWindow[1], "| y =", cutWindow[2])
    if(length(cutWindow) > 3){
      cat(" to a rectangle of a x.width =", cutWindow[3], "and y.height =", cutWindow[4],"\n")
    } else {
      cat(" to a square of a width =", cutWindow[3], "\n")
    }
  }
  
  if(circleRadius != 0){
    cat("Circle radius is", circleRadius, "m\n")
  }
  if(fixedLimit != 0){
    cat(paste0("Plotting only a square from -", fixedLimit, " to +", fixedLimit," m.\n"))
  }
  
  if(quickPreview){
    cat("\n\nONLY FOR QUICK PREVIEW!!!\n\n")
  }
  
  
  
  # functions formzahl and tree species
  {
    treeSpecies <- function(number) {
      species <- "YY"
      
      species <- switch (paste0(number),
                         "1" = "Fi",
                         "2" = "Ta",
                         "3" = "La",
                         "4" = "Ki",
                         "5" = "Sk",
                         "6" = "Zk",
                         "10" = "Bu",
                         "11" = "Ei",
                         "12" = "Ha",
                         "13" = "Es",
                         "14" = "Ah",
                         "15" = "Ul",
                         "18" = "Ks",
                         "19" = "Bi",
                         "20" = "Er",
                         "22" = "Wp",
                         "23" = "Zp",
                         "26" = "We",
                         "28" = "Hs",
                         "NA")
      return(species)
    }
    
    treeSpecies <- Vectorize(treeSpecies)
    
    
    treeSpeciesNumber <- function(species) {
      number <- -1
      
      number <- switch (paste0(species),
                        "Fi" = 1,
                        "Ta" = 2,
                        "La" = 3,
                        "Ki" = 4,
                        "Sk" = 5,
                        "Zk" = 6,
                        "Bu" = 10,
                        "Vb" = 10,
                        "Ei" = 11,
                        "Ze" = 11,
                        "Ha" = 12,
                        "Hb" = 12,
                        "Es" = 13,
                        "Ah" = 14,
                        "Li" = 14,
                        "Ul" = 15,
                        "Ks" = 18,
                        "Bi" = 19,
                        "Se" = 20,
                        "Er" = 20,
                        "Wp" = 22,
                        "Zp" = 23,
                        "Pa" = 23,
                        "We" = 26,
                        "XL" = 26,
                        "Hs" = 28,
                        "Lt" = 28,
                        -1)
      
      return(number)
    }
    
    treeSpeciesNumber <- Vectorize(treeSpeciesNumber)
    
    
    
    
    # Formzahl nach Pollanschuetz
    
    fPoll <- function(species,dbh_cm,height_m) {
      BHD = dbh_cm / 10
      Hoehe = height_m * 10
      
      if (species == 1 && BHD > 1.05) { #Fichte, sonst. Nadelholz
        b1 = 0.46818
        b2 = -0.013919
        b3 = -28.213
        b4 = 0.37474
        b5 = -0.28875
        b6 = 28.279
        b7 = 0
      }
      if (species == 1 && BHD <= 1.05) { #Fichte, sonst. Nadelholz
        b1 = 0.563443
        b2 = -0.12731
        b3 = -8.55022
        b4 = 0
        b5 = 0
        b6 = 7.6331
        b7 = 0
      }
      if (species == 2 && BHD > 1.05) { #Tanne
        b1 = 0.580223
        b2 = -0.0307373
        b3 = -17.1507
        b4 = 0.089869
        b5 = -0.080557
        b6 = 19.661
        b7 = -2.45844
      }
      if (species == 2 && BHD <= 1.05) { #Tanne
        b1 = 0.560673
        b2 = 0.15468
        b3 = -0.65583
        b4 = 0.03321
        b5 = 0
        b6 = 0
        b7 = -0
      }
      if (species == 3 && BHD > 1.05) { #L?rche
        b1 = 0.609443
        b2 = -0.0455748
        b3 = -18.6631
        b4 = -0.248736
        b5 = 0.126594
        b6 = 36.9783
        b7 = -14.204
      }
      if (species == 3 && BHD <= 1.05) { #L?rche
        b1 = 0.48727
        b2 = 0
        b3 = -2.04291
        b4 = 0
        b5 = 0
        b6 = 5.9995
        b7 = 0
      }
      if (species == 4) { #Kiefer
        b1 = 0.435949
        b2 = -0.0149083
        b3 = 5.21091
        b4 = 0
        b5 = 0.028702
        b6 = 0
        b7 = 0
      }
      if (species == 5) { #Schwarzkiefer
        b1 = 0.53438
        b2 = -0.00763
        b3 = 0
        b4 = 0
        b5 = 0
        b6 = 0
        b7 = 2.2414
      }
      if (species == 6) { #Zirbe
        b1 = 0.525744
        b2 = -0.0334896
        b3 = 7.38943
        b4 = -0.10646
        b5 = 0
        b6 = 0
        b7 = 3.34479
      }
      if (species == 10 && BHD > 1.05) { #Buche, Kastanie, Robinie, S|bus
        b1 = 0.686253
        b2 = -0.0371508
        b3 = -31.0674
        b4 = -0.386321
        b5 = 0.219462
        b6 = 49.6163
        b7 = -22.3719
      }
      if (species == 10  && BHD <= 1.05) { #Buche, Kastanie, Robinie, S|bus
        b1 = 0.5173
        b2 = 0
        b3 = -13.62144
        b4 = 0
        b5 = 0
        b6 = 9.9888
        b7 = 0
      }
      if (species == 11 && BHD > 1.05) { #Eiche
        b1 = 0.115631
        b2 = 0
        b3 = 65.9961
        b4 = 1.20321
        b5 = -0.930406
        b6 = -215.758
        b7 = 168.477
      }
      if (species == 11 && BHD <= 1.05) { #Eiche
        b1 = 0.417118
        b2 = 0.21941
        b3 = 13.32594
        b4 = 0
        b5 = 0
        b6 = 0
        b7 = 0
      }
      if (species == 12) { #Hainbuche
        b1 = 0.32473
        b2 = 0.02432
        b3 = 0
        b4 = 0.23972
        b5 = 0
        b6 = -9.9388
        b7 = 0
      }
      if (species == 13) { #Esche
        b1 = 0.48122
        b2 = -0.01489
        b3 = -10.83056
        b4 = 0
        b5 = 0
        b6 = 9.3936
        b7 = 0
      }
      if (species == 14 ) { #Ah|n, Linde
        b1 = 0.50101
        b2 = -0.03521
        b3 = -8.07176
        b4 = 0
        b5 = 0.03521
        b6 = 0
        b7 = 0
      }
      if (species == 15) { #Ulme
        b1 = 0.44215
        b2 = -0.02446
        b3 = 0
        b4 = 0
        b5 = 0
        b6 = 0
        b7 = 2.87714
      }
      if (species == 18) { #Kirsche
        b1 = 0.54008
        b2 = -0.02716
        b3 = -25.11447
        b4 = 0.08327
        b5 = 0
        b6 = 9.3988
        b7 = 0
      }
      if (species == 19) { #Birke
        b1 = 0.42831
        b2 = -0.06643
        b3 = 0
        b4 = 0
        b5 = 0
        b6 = 8.4307
        b7 = 0
      }
      if (species == 20 && BHD > 1.05) { #Erle
        b1 = 0.42937
        b2 = 0
        b3 = -4.10259
        b4 = 0
        b5 = 0
        b6 = 16.7578
        b7 = -5.16631
      }
      if (species == 20 && BHD <= 1.05) { #Erle
        b1 = 0.387399
        b2 = 0
        b3 = 7.17123
        b4 = 0.04407
        b5 = 0
        b6 = 0
        b7 = 0
      }
      if (species == 22 && BHD > 1.05) { #Wei?pappel
        b1 = 0.31525
        b2 = 0
        b3 = 0
        b4 = 0.51079
        b5 = -0.34279
        b6 = -26.08
        b7 = 28.6334
      }
      if (species == 22 && BHD <= 1.05) { #Wei?pappel
        b1 = 0.366419
        b2 = 0
        b3 = 1.13323
        b4 = 0.1306
        b5 = 0
        b6 = 0
        b7 = 0
      }
      if (species == 23 ) {  #Schwarzpappel, Zitterppappel
        b1 = 0.4115
        b2 = -0.00989
        b3 = -28.27478
        b4 = 0.35599
        b5 = -0.21986
        b6 = 21.4913
        b7 = 0
      }
      if (species == 26 ) { #Salix, Kirsche, Other broadleaved trees
        b1 = 0.54008
        b2 = -0.02716
        b3 = -25.11447
        b4 = 0.08327
        b5 = 0
        b6 = 9.3988
        b7 = 0
      }
      if (species == 28 ) { #Hasel, Latsche
        b1 = 0
        b2 = 0
        b3 = 0
        b4 = 0
        b5 = 0
        b6 = 0
        b7 = 0
      }
      
      f = b1 + b2*(log(BHD))^2 + b3*(1/Hoehe) + b4*(1/BHD) + b5*(1/(BHD^2)) + b6*(1/(BHD*Hoehe)) + b7*(1/(BHD^2*Hoehe))
      
      f
      
    }
    
    fPoll <- Vectorize(fPoll)
  }
  
  
  
  appPath <- paste0(dirPath, "app/")
  
  if(!dir.exists(appPath)) dir.create(appPath)
  picPath.app.temp <- paste0(appPath, "temp/")
  if(createBGR_pic){
    if(!dir.exists(picPath.app.temp)) dir.create(picPath.app.temp, recursive = T)
    }

  
  
  if(exists("borderLines")){
    # Convex hull + buffer 0 m radius extend dilation
    # set alpha value quite high, so that the points are spanned and no holes are in between
    borderHull <- ashape(borderLines$x, borderLines$y, alpha = 150)
    owin.window <- owin(xrange=range(borderLines$x), yrange=range(borderLines$y))
    trees_edges_ppp <- psp(borderHull$edges[,3], borderHull$edges[,4],
                           borderHull$edges[,5], borderHull$edges[,6], window=owin.window)
    hull.trees.ext <- dilation(trees_edges_ppp, r=0.00001)
    #plot(hull.trees.ext, border= 2)
    trees_edges_buffer <- owin(poly=hull.trees.ext$bdry[[1]])
    area_m <- round(area.owin(trees_edges_buffer),2)
    cat("Framed area has in total: ", area_m, "m?\n")
  }
  
  
  
  
  
  timeRead1 <- Sys.time()
  # CREATING PNG BGR PICTURE ####
  if(createBGR_pic){
    
    if(!is.na(fileFinder[1])){
      ## a) FROM FILEFINDER ####
      for(i in 1:length(slices)){
        
        #cat("reading in file", paste0(groundPath, fileFinder, slices[i]), "\n")
        cat(paste0("Reading in ", slices[i], "... "))
        for(j in 1:length(fileFinder)){
          if(j == 1){
            co <- capture.output(slice <- readLAS(paste0(groundPath, fileFinder[j], slices[i]), 
                             select = selector))
          } else {
            co <- capture.output(newSlice <- readLAS(paste0(groundPath, fileFinder[j], slices[i]), 
                                                     select = selector))
            slice <- rbind(slice, newSlice)
          }
        }
        if(!colorVersion){
          if(drawRedSlice){
            cat(paste0("_clusterSlice_120to140.laz... "))
            #cat("reading in redSlice", paste0(groundPath, fileFinder, "_clusterSlice_120to140.laz"), "\n")
            for(j in 1:length(fileFinder)){
              if(j == 1){
                co <- capture.output(redSlice <- 
                                       readLAS(paste0(groundPath, 
                                                      fileFinder[j], "_clusterSlice_120to140.laz"), 
                                               select = selector))
              } else {
                co <- capture.output(newRedSlice <- readLAS(paste0(groundPath, fileFinder[j], "_clusterSlice_120to140.laz"), 
                                                             select = selector))
                redSlice <- rbind(redSlice, newRedSlice)
              }
            }
          }
          if(drawGround){
            cat(paste0(slices_gr, "... "))
            #cat("reading in ground", paste0(groundPath, fileFinder, slices_gr), "\n")
            for(j in 1:length(fileFinder)){
              if(j == 1){
                ground <- readLAS(paste0(groundPath, fileFinder[j], slices_gr), 
                                  select = selector)
              } else {
                co <- capture.output(newGroundSlice <- readLAS(paste0(groundPath, fileFinder[j], slices_gr), 
                                                               select = selector))
                ground <- rbind(ground, newGroundSlice)
              }
            }
          }
          if(addTrajLAZ){
            cat(paste0(slices_grTRAJ, "... "))
            #cat("reading in ground", paste0(groundPath, fileFinder, slices_gr), "\n")
            for(j in 1:length(fileFinder)){
              if(j == 1){
                trajLAZ <- readLAS(paste0(groundPath, fileFinder[j], slices_grTRAJ), 
                                   select = selector)
              } else {
                co <- capture.output(newTrajLAZ <- readLAS(paste0(groundPath, fileFinder[j], slices_grTRAJ), 
                                      select = selector))
                trajLAZ <- rbind(trajLAZ, newTrajLAZ)
              }
            }
          }
          if (sum(cutWindow == c(-1000, -1000, 2000, 2000)) != 4) {
            cat("\nPre-cutting slice ")
            
            XL <- cutWindow[1]
            YL <- cutWindow[2]
            width <- cutWindow[3]
            if(length(cutWindow)>3){
              height <- cutWindow[4]
            } else {
              height <- width
            }
            slice <- filter_poi(slice, X > XL, X < XL + width, 
                                Y > YL, Y < YL + height)
            if(drawRedSlice){
              cat("- red slice ")
              redSlice <- filter_poi(redSlice, X > XL, X < XL + width, 
                                     Y > YL, Y < YL + height)
            } 
            if(drawGround){
              cat("- ground ")
              ground <- filter_poi(ground, X > XL, X < XL + width, 
                                   Y > YL, Y < YL + height)
            } 
            if(addTrajLAZ){
              cat("- trajlaz ")
              trajLAZ <- filter_poi(trajLAZ, X > XL, X < XL + width, 
                                    Y > YL, Y < YL + height)
            }
            cat("done!\n")
          }
          slice_sm <- decimate_points(slice, random(100))
          
        }
        timeRead2 <- Sys.time()
        cat("done.\n")
        cat("All files read in a ")
        print.difftime(round(timeRead2 - timeRead1), 1)
        
        
        
        
        cat("\nNow is", format(Sys.time()), "\n")
        ### dimensioning ####
        # get dimensions of the first slice (span) rounded to 10 meters (30 m or 40 m or 1670 m :D)
        if(i == 1){
          minY <- slice@header@PHB$`Min Y`
          maxY <- slice@header@PHB$`Max Y`
          minX <- slice@header@PHB$`Min X`
          maxX <- slice@header@PHB$`Max X`
          
          minY <- floor(minY/10)*10
          minX <- floor(minX/10)*10
          maxY <- ceiling(maxY/10)*10
          maxX <- ceiling(maxX/10)*10
          
          if(fixedLimit != 0){
            minY <- -fixedLimit
            minX <- -fixedLimit
            maxY <- fixedLimit
            maxX <- fixedLimit
          }
          
          spanY <- maxY - minY
          spanX <- maxX - minX
          spanY_cm <- spanY * 100
          spanX_cm <- spanX * 100
          
          
          topLeftY <- ceiling(maxY * 100)
          topLeftX <- floor(minX * 100)
          bottomRightY <- floor(minY * 100)
          bottomRightX <- ceiling(maxX * 100)
          
          spanY_cm_drawn <- topLeftY - bottomRightY
          spanX_cm_drawn <- bottomRightX - topLeftX
          
          
          # RESOLUTION PREVIEW
          if(1 == 2){
            {
              ### resolution preview
              prevXLeft <- 0
              prevXRight <- 5
              prevYDown <- 0
              prevYUp <- 5
              slicePrev <- filter_poi(slice, X < prevXRight, X > prevXLeft, Y > prevYDown, Y < prevYUp)
              #plot(slicePrev)
              png(paste0(picPath.app.temp, fileFinder, slices[i],"_prev_",pixelUnit_cm,"px_cm.png"),
                  width = (prevXRight - prevXLeft)*100*pixelUnit_cm, height = (prevYUp - prevYDown)*100*pixelUnit_cm, type = "cairo")
              par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
              plot(1, type = "n", xlim = c(prevXLeft, prevXRight), ylim = c(prevYDown, prevYUp),
                   asp = 1, axes = T)
              maxInt <- max(slicePrev@data$Intensity)
              slicePrev@data$col <- as.integer((slicePrev@data$Intensity/maxInt) * 255)
              range(slicePrev@data$col)
              points(slicePrev@data$X, slicePrev@data$Y, cex = 0.00001,
                     col = rgb(0,0,0, alpha = slicePrev@data$Intensity, maxColorValue = maxInt))
              points(slicePrev@data$X, slicePrev@data$Y, cex = 0.00001,
                     col = rgb(0,0,0, alpha = slicePrev@data$Intensity, maxColorValue = maxInt))
              
              
              dev.off()
            }
          }
        }
        
        
        #spanX_cm_drawn <- 42000
        #spanY_cm_drawn <- 30000
        #pixelUnit_cm <- 1
        xLarge <- (spanX_cm_drawn*pixelUnit_cm) > 32768
        yLarge <- (spanY_cm_drawn*pixelUnit_cm) > 32768
        areaLarge <- (spanY_cm_drawn*pixelUnit_cm)*(spanX_cm_drawn*pixelUnit_cm) > 1000000000
        
        while(xLarge | yLarge | areaLarge){
          cat("Cutting pixelUnit to", pixelUnit_cm - 0.05, "because of")
          if(xLarge) cat(" x =", (spanX_cm_drawn*pixelUnit_cm))
          if(yLarge) cat(" y =", (spanY_cm_drawn*pixelUnit_cm))
          if(areaLarge) cat(" area =", round((spanY_cm_drawn*pixelUnit_cm)*(spanX_cm_drawn*pixelUnit_cm)/1000000), "Mio")
          
          pixelUnit_cm <- pixelUnit_cm - 0.05
          xLarge <- (spanX_cm_drawn*pixelUnit_cm) > 32768
          yLarge <- (spanY_cm_drawn*pixelUnit_cm) > 32768
          areaLarge <- (spanY_cm_drawn*pixelUnit_cm)*(spanX_cm_drawn*pixelUnit_cm) > 1000000000
          cat("\n")
        }
        
        
        
        if(pixelUnit_cm < 0.5){
          cat("Reducing all input files - ")
          if(drawGround){
            cat("ground ")
            tGround <- ground@header@PHB$`Number of point records`
            ground <- decimate_points(ground, random(400))
            tGround <- ground@header@PHB$`Number of point records` / tGround
            cat("x", round(tGround,2), " - ", sep = "")
          }
          if(drawRedSlice){
            cat("red ")
            tRed <- redSlice@header@PHB$`Number of point records`
            redSlice <- decimate_points(redSlice, random(300))
            tRed <- redSlice@header@PHB$`Number of point records` / tRed
            cat("x", round(tRed,2), " - ", sep = "")
          }
          cat("black ")
          tSlice <- slice@header@PHB$`Number of point records` 
          slice <- decimate_points(slice, random(1000))
          tSlice <- slice@header@PHB$`Number of point records` / tSlice
          cat("x", round(tSlice,2), " - ", sep = "")
          cat("done!\n")
        }
        
        
        
        ## plotting background image of fileFinders ####
        {
          bgrTime <- Sys.time()
          cat("EXHAUSTIVE: plot background image",spanX_cm_drawn*pixelUnit_cm,"x",spanY_cm_drawn*pixelUnit_cm,"pixel (area =", round(spanX_cm_drawn*pixelUnit_cm * spanY_cm_drawn*pixelUnit_cm / 1000000), "Mio)...\n")
          cat("            resolution is", pixelUnit_cm, "pixel per cm - using slice", gsub(".laz", "", slices[i]), "\n")
          cat("            center point 0 | 0 will be drawn at: \n")
          cat(rep(" ", 26 - nchar(-topLeftX*pixelUnit_cm)), (-topLeftX)*pixelUnit_cm," | ",topLeftY*pixelUnit_cm," pixel count in the image.\n", sep = "")
          cat("(takes time...) ")
        }
        tryCatch(
          {
            #png(paste0(picPath.app.temp, fileFinder, "_", clip.radius,slices[i],".png"), width = 8000, height = 8000, type = "cairo")
            
            if(createTIFF){
              tiff(paste0(picPath.app.temp, outName, slices[i],".png"), 
                  width = spanX_cm_drawn*pixelUnit_cm, height = spanY_cm_drawn*pixelUnit_cm, type = "cairo")
              
            } else {
              png(paste0(picPath.app.temp, outName, slices[i],".png"),
                  width = spanX_cm_drawn*pixelUnit_cm, height = spanY_cm_drawn*pixelUnit_cm, type = "cairo")
            }
            par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
            plot(1, type = "n", xlim = c(topLeftX/100,bottomRightX/100), ylim = c(bottomRightY/100, topLeftY/100),
                 asp = 1, axes = T)
            
            
            
            
            
            if(colorVersion){
              
              
              #better, but didnt work for JA7
              #groundAllPath <- paste0(groundPath, fileFinder, "_clusterSlice_50to200.laz")
              groundAllPath <- paste0(groundPath, fileFinder, "_ground.laz")
              cat("Reading in complete ground", groundAllPath, "\n")
              co <- capture.output(groundAll <- readLAS(groundAllPath, select = selector))
              cat("Merging ground and slice...")
              bothLAS <- rbind(groundAll, slice)
              cat(" and filter duplicates...\n")
              bothLAS <- filter_duplicates(bothLAS)
              # NORMALIZE HEIGHT
              cat(" - done, reading dtm ")
              dtm <- raster(groundModel.path)
              
              
              if (sum(cutWindow == c(-1000, -1000, 2000, 2000)) != 4) {
                cat("Cutting dtm... ")
                dtmPuffer <- 2 #m for height lines to not get split
                if(length(cutWindow) > 3){
                  e <- c(cutWindow[1]-dtmPuffer, 
                         cutWindow[2]-dtmPuffer, 
                         cutWindow[1]+cutWindow[3]+dtmPuffer, 
                         cutWindow[2]+cutWindow[4]+dtmPuffer)
                } else {
                  e <- c(cutWindow[1]-dtmPuffer,
                         cutWindow[2]-dtmPuffer, 
                         cutWindow[1]+cutWindow[3]+dtmPuffer, 
                         cutWindow[2]+cutWindow[3]+dtmPuffer)
                }
                
                dtm <- crop(dtm, e)
              }
              #crs(dtm) <- CRS("+init=epsg:31258")
              cat("to normalize height...")
              bothLAS <- normalize_height(bothLAS, dtm, na.rm = TRUE)
              cat("ok.\n")
              
              # To-Do: if all slices are present, then dont do normalisation (takes a lot of time)
              
              
              slicePath.red <- paste0(groundPath, fileFinder, "_xSlice_BHD_red.laz")
              if(file.exists(slicePath.red)){
                co <- capture.output(sliceRed <- readLAS(slicePath.red, select = selector))
                cat("Red slice read.\n")
              } else {
                if(thin){
                  sliceRed <- filter_poi(bothLAS, Z < 1.33, Z > 1.27)
                } else {
                  sliceRed <- filter_poi(bothLAS, Z < 1.4, Z > 1.2)
                }
                sliceRed <- unnormalize_height(sliceRed)
                writeLAS(sliceRed, slicePath.red)
                cat("Red slice created.\n")
              }
              
              
              
              slicePath.brown <- paste0(groundPath, fileFinder, "_xSlice_100to110_brown.laz")
              if(file.exists(slicePath.brown)){
                co <- capture.output(sliceBrown <- readLAS(slicePath.brown, select = selector))
                cat("Brown slice read.\n")
              } else {
                sliceBrown <- filter_poi(bothLAS, Z <= 1.1, Z > 1.0)
                sliceBrown <- unnormalize_height(sliceBrown)
                writeLAS(sliceBrown, slicePath.brown)
                cat("Brown slice created.\n")
              }
              
              slicePath.blue <- paste0(groundPath, fileFinder, "_xSlice_190to200_blue.laz")
              if(file.exists(slicePath.blue)){
                co <- capture.output(sliceBlue <- readLAS(slicePath.blue, select = selector))
                cat("Blue slice read.\n")
              } else {
                sliceBlue <- filter_poi(bothLAS, Z <= 2, Z > 1.9)
                sliceBlue <- unnormalize_height(sliceBlue)
                writeLAS(sliceBlue, slicePath.blue)
                cat("Blue slice created.\n")
              }
              
              slicePath.grey <- paste0(groundPath, fileFinder, "_xSlice_200to300_grey.laz")
              if(file.exists(slicePath.grey)){
                co <- capture.output(sliceGrey <- readLAS(slicePath.grey, select = selector))
                cat("Grey slice read.\n")
              } else {
                sliceGrey <- filter_poi(bothLAS, Z > 2, Z <= 3)
                sliceGrey <- unnormalize_height(sliceGrey)
                writeLAS(sliceGrey, slicePath.grey)
                cat("Grey slice created.\n")
              }
              
              slicePath.black <- paste0(groundPath, fileFinder, "_xSlice_110to190_black.laz")
              if(file.exists(slicePath.black)){
                co <- capture.output(sliceBlack <- readLAS(slicePath.black, select = selector))
                cat("Black slice read.\n")
              } else {
                sliceBlack <- filter_poi(bothLAS, Z > 1.1, Z <= 1.9)
                sliceBlack <- unnormalize_height(sliceBlack)
                writeLAS(sliceBlack, slicePath.black)
                cat("Black slice created.\n")
              }
              
              slicePath.green <- paste0(groundPath, fileFinder, "_xSlice_050to100_green.laz")
              if(file.exists(slicePath.green)){
                co <- capture.output(sliceGreen <- readLAS(slicePath.green, select = selector))
                cat("Green slice read.\n")
              } else {
                sliceGreen <- filter_poi(bothLAS, Z > 0.5, Z <= 1)
                sliceGreen <- unnormalize_height(sliceGreen)
                writeLAS(sliceGreen, slicePath.green)
                cat("Green slice created.\n")
              }
              
              
              cat("\nMixing on our color palette: \n")
              maxInt <- max(sliceGreen@data$Intensity)
              points(sliceGreen@data$X, sliceGreen@data$Y, cex = 0.00001,
                     col = rgb(0,maxInt,0, alpha = sliceGreen@data$Intensity, maxColorValue = maxInt))
              cat("green")
              
              maxInt <- max(sliceGrey@data$Intensity)
              points(sliceGrey@data$X, sliceGrey@data$Y, cex = 0.00001,
                     col = rgb(maxInt/2, maxInt/2, maxInt/2, alpha = sliceGrey@data$Intensity, maxColorValue = maxInt))
              cat(" - grey")
              
              maxInt <- max(sliceBlack@data$Intensity)/1.5
              sliceBlack@data$Intensity[sliceBlack@data$Intensity > maxInt] <- as.integer(maxInt)
              points(sliceBlack@data$X, sliceBlack@data$Y, cex = 0.00001,
                     col = rgb(0,0,0, alpha = sliceBlack@data$Intensity, maxColorValue = maxInt))
              cat(" - black")
              
              
              maxInt <- max(sliceBrown@data$Intensity)
              points(sliceBrown@data$X, sliceBrown@data$Y, cex = 0.00001,
                     col = rgb(maxInt/3,maxInt/4,0, alpha = sliceBrown@data$Intensity, maxColorValue = maxInt))
              cat(" - brown")
              
              
              maxInt <- max(sliceBlue@data$Intensity)
              points(sliceBlue@data$X, sliceBlue@data$Y, cex = 0.00001,
                     col = rgb(0, 0, maxInt, alpha = sliceBlue@data$Intensity, maxColorValue = maxInt))
              cat(" - blue")
              
              maxInt <- max(sliceRed@data$Intensity)
              points(sliceRed@data$X, sliceRed@data$Y, cex = 0.00001,
                     col = rgb(maxInt, 0, 0, alpha = sliceRed@data$Intensity, maxColorValue = maxInt))
              cat(" - and red!\n")
              
              
            } else {
              # Traditionally with only the red band at DBH
              # no nice colors
              
              maxInt <- as.integer(max(slice@data$Intensity)/1.5)
              slice@data$Intensity[slice@data$Intensity > maxInt] <- as.integer(maxInt)
              
              cat(".")
              if(drawGround){
                maxIntGr <- as.integer(max(ground@data$Intensity)/1.5)
                ground@data$Intensity[ground@data$Intensity > maxIntGr] <- as.integer(maxIntGr)
                if(fullGroundGrey){
                  cat("0")
                  points(ground@data$X, ground@data$Y, cex = 0.00001, 
                         col = rgb(maxIntGr/2, maxIntGr/2, maxIntGr/2, 
                                   alpha = ground@data$Intensity, maxColorValue = maxIntGr))
                } else {
                  cat("o")
                  points(ground@data$X, ground@data$Y, cex = 0.00001, 
                         col = rgb(0,0,0, alpha = ground@data$Intensity, maxColorValue = maxIntGr))
                }
              }
              
              if(addTrajLAZ){
                cat("+TRAJLAZ")
                maxIntTraj <- max(trajLAZ@data$Intensity, na.rm = T)
                plot(trajLAZ@data$X, trajLAZ@data$Y,  cex = pixelUnit_cm, #pch = 16,
                     col = rgb(0, 0, maxIntTraj, alpha = trajLAZ@data$Intensity, maxColorValue = maxIntTraj))
              }
              
              if(drawTraj){
                cat("+TRAJ")
                points(allTraj$y ~ allTraj$x, col = allTraj$col, cex = pixelUnit_cm, pch = 16)
              }
              if(nrow(greySpots) > 0){
                cat(paste0("+",nrow(greySpots),"spots"))
                if(ncol(greySpots) > 2){
                  if(ncol(greySpots) == 3){
                    greySpots[,4] <- rgb(200,200,200, alpha = 100, maxColorValue = 255)
                  }
                  points(greySpots[,1], greySpots[,2], pch = 19, 
                         cex = greySpots[,3]*pixelUnit_cm*20, col = greySpots[,4])
                  
                }
                
              }
              if(nrow(drawLines) > 0){
                cat(paste0("+",nrow(drawLines),"lines"))
                if(ncol(drawLines) > 4){
                  if(ncol(drawLines) == 4){
                    drawLines[,5] <- rgb(200,200,200, alpha = 100, maxColorValue = 255)
                  }
                  if(ncol(drawLines) == 5){
                    drawLines[,6] <- 1 # m thick lines
                  }
                  segments(x0 = drawLines[,1], y0 = drawLines[,2], 
                           x1 = drawLines[,3], y1 = drawLines[,4], 
                           col = drawLines[,5], 
                           lwd = drawLines[,6]*pixelUnit_cm/3*400)
                  
                  
                }
                
              }
              
              
              if(isoLines > 0){
                cat("+ISO-LINES", isoLines, "m")
                if(file.exists(dtm.path)){
                  dtm <- raster(dtm.path)
                  if (sum(cutWindow == c(-1000, -1000, 2000, 2000)) != 4) {
                    cat(" CROPPED")
                    dtmPuffer <- 2 #m for height lines to not get split
                    if(length(cutWindow) > 3){
                      e <- c(cutWindow[1]-dtmPuffer, 
                             cutWindow[2]-dtmPuffer, 
                             cutWindow[1]+cutWindow[3]+dtmPuffer, 
                             cutWindow[2]+cutWindow[4]+dtmPuffer)
                    } else {
                      e <- c(cutWindow[1]-dtmPuffer,
                             cutWindow[2]-dtmPuffer, 
                             cutWindow[1]+cutWindow[3]+dtmPuffer, 
                             cutWindow[2]+cutWindow[3]+dtmPuffer)
                    }
                  }
                  dtm <- crop(dtm, e)
                  #crs(dtm) <- CRS("+init=epsg:31258")
                  
                  dtm <- disaggregate(dtm, 2, method='bilinear')
                  dtm <- focal(dtm, w = focalWeight(dtm, 0.9, "Gauss"))
                  
                  intervalLines <- isoLines
                  linesInt <- seq(floor(dtm@data@min/intervalLines)*intervalLines, 
                                  ceiling(dtm@data@max/intervalLines)*intervalLines, by = intervalLines)
                  
                  intervalLines5 <- isoLines * 5
                  linesInt5 <- seq(floor(dtm@data@min/intervalLines5)*intervalLines5, 
                                   ceiling(dtm@data@max/intervalLines5)*intervalLines5, by = intervalLines5)
                  
                  options("max.contour.segments" = 100000)
                  
                  
                  useLWD <- 1.6
                  if(pixelUnit_cm <= 0.8) useLWD <- 2.5
                  #maxCells gives us the area of the plot, so that we have 
                  # one edge of the contours every 10 cm
                  maxCells <- round(100*(spanY_cm_drawn*pixelUnit_cm)*(spanX_cm_drawn*pixelUnit_cm))
                  terra::contour(dtm, maxcells = maxCells, levels = linesInt, 
                                 add=TRUE, col = "orange", labcex = 0.0001, lwd = useLWD)
                  terra::contour(dtm, maxcells = maxCells, levels = linesInt5, 
                                 add=TRUE, col = "orange", labcex = 0.0001, lwd = useLWD + 4.5)
                  
                  cat(" in orange+")
                } else {
                  cat("NO DTM MODEL SPECIFIED - SKIPPING... HUOM - Implement it from fileFinder!")
                }
                
              } else {
                cat(".")
              }
              
              
              
              if(quickPreview){
                points(slice_sm@data$X, slice_sm@data$Y, cex = 0.00001, col = rgb(0,0,0, alpha = slice_sm@data$Intensity, maxColorValue = maxInt))
              } else {
                points(slice@data$X, slice@data$Y, cex = 0.00001, col = rgb(0,0,0, alpha = slice@data$Intensity, maxColorValue = maxInt))
              }
              
              cat(".")
              if(drawRedSlice){
                maxInt <- max(redSlice@data$Intensity)
                points(redSlice@data$X, redSlice@data$Y, cex = 0.00001, col = rgb(maxInt,0,0, alpha = redSlice@data$Intensity, maxColorValue = maxInt))
              }
            }
            
            
            
            
            
            #draw.circle(0,0,40)
            
            
            
            # PLOT ROSALIA RIVONA BORDERS
            if(circleRadius == 0 && exists("borderLines")){
              #corners <- data.frame("x" = c(1.13413, -3.56123, 157.593, 111.246, 1.13413), "y" = c(-0.0310049, 76.4206, -16.573, -63.6488, -0.0310049))
              corners <- data.frame("x" = c(borderLines$x,borderLines$x[1]) , "y" = c(borderLines$y,borderLines$y[1]))
              points(corners$x, corners$y, pch = 16, cex = 5, col = "red")
              lines(corners$x, corners$y, lty = 2, lwd = 2, col = "red")
            }
            
            if (greenTreeLocations) {
              treeListFile <- paste0(dirPath,fileFinder,"/trees_dbh.txt")
              if(!file.exists(treeListFile[1])){
                treeListFile <- paste0(dirPath,fileFinder,"_ALLGO_100to300_dbh/trees_dbh.txt")
              }
              for(j in 1:length(fileFinder)){
                if(file.exists(treeListFile[i])){
                  trees <- read.table(treeListFile[i], header = TRUE)
                  #trees$dCol <- trees$dbh/max(trees$dbh)
                  trees$dCol <- rgb(50,250,50, maxColorValue = 255)
                  trees$dCol[trees$id >= 9000] <- rgb(255,200,0, maxColorValue = 255)
                  
                  dbh.scaler <- 5
                  
                  points(trees$x, trees$y, pch = 16, cex = log(trees$dbh)/log(dbh.scaler)/1.5, col = trees$dCol)
                  
                }
              }
            }
            
            if(circleRadius > 0){
              draw.circle(0, 0, radius = circleRadius, lwd = 5, border = "blue")
            }
            
            # central cross
            clip(x1 = -1, x2 = 1,
                 y1 = -2, y2 = 3)
            abline(v=0, h=0, lwd = 2)
            
          }, error = function(error_condition) {
            cat("Error in creating this image, next loop!\n")
            next
          })
      
        dev.off()
        }
      
    }
    else {
      ## b) FROM LAS FILE ####
      # reading in many las files, merging to one
      cat("Reading in file ")
      slice <- ""
      for(i in 1:length(laz.path)){
        cat(i, "- ")
        co <- capture.output(tempLAS <- readLAS(laz.path[i], select = selector))
        #tempLAS <- lasaddextrabytes(tempLAS, i, name = "Part Name", desc = "Segmented part number, Ra = 1 south, Rf = 6 north")
        
        if(i == 1){
          slice <- tempLAS
        } else {
          slice <- rbind(slice, tempLAS)
        }
      }
      cat("done!\n")
      cat("Reading DTM model from", basename(groundModel.path))
      dtm <- raster(groundModel.path)
      cat(" - done, normalizing height...")
      slice <- normalize_height(slice, dtm, na.rm = TRUE)
      cat("ok.\n")
      slice_sm <- decimate_points(slice, random(100))
      if(thin){
        redSlice <- filter_poi(slice, Z < 1.33, Z > 1.27)
      } else {
        redSlice <- filter_poi(slice, Z < 1.4, Z > 1.2)
      }
      redSlice <- unnormalize_height(redSlice)
      
      if(!is.na(slice.cut.low)){
        slice <- filter_poi(slice, Z <= slice.cut.high, Z >= slice.cut.low)
      }
      
      
      ### dimensioning ####
      {
        # get dimensions of the first slice (span) rounded to 10 meters (30 m or 40 m or 1670 m :D)
        minY <- slice@header@PHB$`Min Y`
        maxY <- slice@header@PHB$`Max Y`
        minX <- slice@header@PHB$`Min X`
        maxX <- slice@header@PHB$`Max X`
        
        minY <- floor(minY/10)*10
        minX <- floor(minX/10)*10
        maxY <- ceiling(maxY/10)*10
        maxX <- ceiling(maxX/10)*10
        
        if(fixedLimit != 0){
          minY <- -fixedLimit
          minX <- -fixedLimit
          maxY <- fixedLimit
          maxX <- fixedLimit
        }
        # quick fix to make it 6000 x 6000 px for standard output files
        
        spanY <- maxY - minY
        spanX <- maxX - minX
        spanY_cm <- spanY * 100
        spanX_cm <- spanX * 100
        
        
        topLeftY <- ceiling(maxY * 100)
        topLeftX <- floor(minX * 100)
        bottomRightY <- floor(minY * 100)
        bottomRightX <- ceiling(maxX * 100)
        
        spanY_cm_drawn <- topLeftY - bottomRightY
        spanX_cm_drawn <- bottomRightX - topLeftX
        
        
        # RESOLUTION PREVIEW
        if(1 == 2){
          {
            ### resolution preview
            prevXLeft <- 0
            prevXRight <- 5
            prevYDown <- 0
            prevYUp <- 5
            slicePrev <- filter_poi(slice, X < prevXRight, X > prevXLeft, Y > prevYDown, Y < prevYUp)
            #plot(slicePrev)
            png(paste0(picPath.app.temp, basename(laz.path), slices[1],"_prev_",pixelUnit_cm,"px_cm.png"),
                width = (prevXRight - prevXLeft)*100*pixelUnit_cm, height = (prevYUp - prevYDown)*100*pixelUnit_cm, type = "cairo")
            par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
            plot(1, type = "n", xlim = c(prevXLeft, prevXRight), ylim = c(prevYDown, prevYUp),
                 asp = 1, axes = T)
            maxInt <- max(slicePrev@data$Intensity)
            slicePrev@data$col <- as.integer((slicePrev@data$Intensity/maxInt) * 255)
            range(slicePrev@data$col)
            points(slicePrev@data$X, slicePrev@data$Y, cex = 0.00001,
                   col = rgb(0,0,0, alpha = slicePrev@data$Intensity, maxColorValue = maxInt))
            points(slicePrev@data$X, slicePrev@data$Y, cex = 0.00001,
                   col = rgb(0,0,0, alpha = slicePrev@data$Intensity, maxColorValue = maxInt))
            
            
            dev.off()
          }
        }
        
      }
      
      bgrTime <- Sys.time()
      for(i in 1:length(slices)){
        ### plotting background image of merged las ####
        cat("EXHAUSTIVE: plot background image",spanX_cm_drawn*pixelUnit_cm,"x",spanY_cm_drawn*pixelUnit_cm,"pixel (area =", round(spanX_cm_drawn*pixelUnit_cm * spanY_cm_drawn*pixelUnit_cm / 1000000), "Mio)...\n")
        cat("            resolution is", pixelUnit_cm, "pixel per cm - using slice", gsub(".laz", "", slices[i]), "\n")
        cat("            center point 0 | 0 will be drawn at: \n")
        cat(rep(" ", 26 - nchar(-topLeftX*pixelUnit_cm)), (-topLeftX)*pixelUnit_cm," | ",topLeftY*pixelUnit_cm," pixel count in the image.\n", sep = "")
        tryCatch(
          {
            #png(paste0(picPath.app.temp, fileFinder, "_", clip.radius,slices[i],".png"), width = 8000, height = 8000, type = "cairo")
            
            
            tryCatch({
              png(paste0(picPath.app.temp, basename(laz.path), slices[i],".png"),
                  width = spanX_cm_drawn*pixelUnit_cm, height = spanY_cm_drawn*pixelUnit_cm, type = "cairo")
              
            }, error = function(error_condition) {
              cat("Error with that big resolution, halving it!\n")
              pixelUnit_cm <- pixelUnit_cm/2
              cat("New pixelUnit is", pixelUnit_cm, "cm\n")
              png(paste0(picPath.app.temp, basename(laz.path), slices[i],".png"),
                  width = spanX_cm_drawn*pixelUnit_cm, height = spanY_cm_drawn*pixelUnit_cm, type = "cairo")
            })
            
            par(mar=c(0, 0, 0, 0), xaxs='i', yaxs='i')
            plot(1, type = "n", xlim = c(topLeftX/100,bottomRightX/100), ylim = c(bottomRightY/100, topLeftY/100),
                 asp = 1, axes = T)
            
            
            if(colorVersion){
              
              
              #better, but didnt work for JA7
              #groundAllPath <- paste0(groundPath, fileFinder, "_clusterSlice_50to200.laz")
              
              cat("In this version we use the input file for the color!\n")
              bothLAS <- (slice)
              
              
              # To-Do: if all slices are present, then dont do normalisation (takes a lot of time)
              
              slicePath.red <- paste0(groundPath, basename(laz.path), "_xSlice_BHD_red.laz")
              if(file.exists(slicePath.red)){
                co <- capture.output(sliceRed <- readLAS(slicePath.red, select = selector))
                cat("Red slice read.\n")
              } else {
                if(thin){
                  sliceRed <- filter_poi(bothLAS, Z < 1.33, Z > 1.27)
                } else {
                  sliceRed <- filter_poi(bothLAS, Z < 1.4, Z > 1.2)
                }
                sliceRed <- unnormalize_height(sliceRed)
                writeLAS(sliceRed, slicePath.red)
                cat("Red slice created.\n")
              }
              
              
              
              slicePath.brown <- paste0(groundPath, basename(laz.path), "_xSlice_100to110_brown.laz")
              if(file.exists(slicePath.brown)){
                co <- capture.output(sliceBrown <- readLAS(slicePath.brown, select = selector))
                cat("Brown slice read.\n")
              } else {
                sliceBrown <- filter_poi(bothLAS, Z <= 1.1, Z > 1.0)
                sliceBrown <- unnormalize_height(sliceBrown)
                writeLAS(sliceBrown, slicePath.brown)
                cat("Brown slice created.\n")
              }
              
              slicePath.blue <- paste0(groundPath, basename(laz.path), "_xSlice_190to200_blue.laz")
              if(file.exists(slicePath.blue)){
                co <- capture.output(sliceBlue <- readLAS(slicePath.blue, select = selector))
                cat("Blue slice read.\n")
              } else {
                sliceBlue <- filter_poi(bothLAS, Z <= 2, Z > 1.9)
                sliceBlue <- unnormalize_height(sliceBlue)
                writeLAS(sliceBlue, slicePath.blue)
                cat("Blue slice created.\n")
              }
              
              slicePath.grey <- paste0(groundPath, basename(laz.path), "_xSlice_200to300_grey.laz")
              if(file.exists(slicePath.grey)){
                co <- capture.output(sliceGrey <- readLAS(slicePath.grey, select = selector))
                cat("Grey slice read.\n")
              } else {
                sliceGrey <- filter_poi(bothLAS, Z > 2, Z <= 3)
                sliceGrey <- unnormalize_height(sliceGrey)
                writeLAS(sliceGrey, slicePath.grey)
                cat("Grey slice created.\n")
              }
              
              slicePath.black <- paste0(groundPath, basename(laz.path), "_xSlice_110to190_black.laz")
              if(file.exists(slicePath.black)){
                co <- capture.output(sliceBlack <- readLAS(slicePath.black, select = selector))
                cat("Black slice read.\n")
              } else {
                sliceBlack <- filter_poi(bothLAS, Z > 1.1, Z <= 1.9)
                sliceBlack <- unnormalize_height(sliceBlack)
                writeLAS(sliceBlack, slicePath.black)
                cat("Black slice created.\n")
              }
              
              slicePath.green <- paste0(groundPath, basename(laz.path), "_xSlice_050to100_green.laz")
              if(file.exists(slicePath.green)){
                co <- capture.output(sliceGreen <- readLAS(slicePath.green, select = selector))
                cat("Green slice read.\n")
              } else {
                sliceGreen <- filter_poi(bothLAS, Z > 0.5, Z <= 1)
                sliceGreen <- unnormalize_height(sliceGreen)
                writeLAS(sliceGreen, slicePath.green)
                cat("Green slice created.\n")
              }
              
              
              cat("\nMixing on our color palette: \n")
              maxInt <- max(sliceGreen@data$Intensity)
              points(sliceGreen@data$X, sliceGreen@data$Y, cex = 0.00001,
                     col = rgb(0,maxInt,0, alpha = sliceGreen@data$Intensity, maxColorValue = maxInt))
              cat("green")
              
              maxInt <- max(sliceGrey@data$Intensity)
              points(sliceGrey@data$X, sliceGrey@data$Y, cex = 0.00001,
                     col = rgb(maxInt/2, maxInt/2, maxInt/2, alpha = sliceGrey@data$Intensity, maxColorValue = maxInt))
              cat(" - grey")
              
              maxInt <- max(sliceBlack@data$Intensity)/1.5
              sliceBlack@data$Intensity[sliceBlack@data$Intensity > maxInt] <- as.integer(maxInt)
              points(sliceBlack@data$X, sliceBlack@data$Y, cex = 0.00001,
                     col = rgb(0,0,0, alpha = sliceBlack@data$Intensity, maxColorValue = maxInt))
              cat(" - black")
              
              
              maxInt <- max(sliceBrown@data$Intensity)
              points(sliceBrown@data$X, sliceBrown@data$Y, cex = 0.00001,
                     col = rgb(maxInt/3,maxInt/4,0, alpha = sliceBrown@data$Intensity, maxColorValue = maxInt))
              cat(" - brown")
              
              
              maxInt <- max(sliceBlue@data$Intensity)
              points(sliceBlue@data$X, sliceBlue@data$Y, cex = 0.00001,
                     col = rgb(0, 0, maxInt, alpha = sliceBlue@data$Intensity, maxColorValue = maxInt))
              cat(" - blue")
              
              maxInt <- max(sliceRed@data$Intensity)
              points(sliceRed@data$X, sliceRed@data$Y, cex = 0.00001,
                     col = rgb(maxInt, 0, 0, alpha = sliceRed@data$Intensity, maxColorValue = maxInt))
              cat(" - and red!\n")
              
              
            } else {
              maxInt <- max(slice@data$Intensity)/1.5
              slice@data$Intensity[slice@data$Intensity > maxInt] <- as.integer(maxInt)
              
              # Traditionally with only the red band at DBH
              if(drawGround){
                maxIntGr <- max(ground@data$Intensity)
                points(ground@data$X, ground@data$Y, cex = 0.00001, col = rgb(0,0,0, alpha = ground@data$Intensity, maxColorValue = maxIntGr))
              }
              if(addTrajLAZ){
                cat("+TRAJLAZ")
                maxIntTraj <- max(trajLAZ@data$Intensity, na.rm = T)
                points(trajLAZ@data$X, trajLAZ@data$Y, cex = pixelUnit_cm, #pch = 16,
                       col = rgb(0, 0, maxIntTraj, alpha = trajLAZ@data$Intensity, maxColorValue = maxIntTraj))
              }
              if(drawTraj){
                cat("+TRAJ")
                points(allTraj$y ~ allTraj$x, col = allTraj$col, cex = pixelUnit_cm, pch = 16)
              }
              if(nrow(greySpots) > 0){
                cat(paste0("+",nrow(greySpots),"spots"))
                if(ncol(greySpots) > 2){
                  if(ncol(greySpots) == 3){
                    greySpots[,4] <- rgb(200,200,200, alpha = 100, maxColorValue = 255)
                  }
                  points(greySpots[,1], greySpots[,2], pch = 19, 
                         cex = greySpots[,3]*pixelUnit_cm*20, col = greySpots[,4])
                  
                }
                
              }
              
              if(quickPreview){
                points(slice_sm@data$X, slice_sm@data$Y, cex = 0.00001, col = rgb(0,0,0, alpha = slice_sm@data$Intensity, maxColorValue = maxInt))
              } else {
                points(slice@data$X, slice@data$Y, cex = 0.00001, col = rgb(0,0,0, alpha = slice@data$Intensity, maxColorValue = maxInt))
              }
              
              if(drawRedSlice){
                maxInt <- max(redSlice@data$Intensity)
                points(redSlice@data$X, redSlice@data$Y, cex = 0.00001, col = rgb(maxInt,0,0, alpha = redSlice@data$Intensity, maxColorValue = maxInt))
              }
            }
            
            
            
            
            
            #draw.circle(0,0,40)
            
            
            
            # PLOT ROSALIA RIVONA BORDERS
            if(circleRadius == 0 && exists("borderLines")){
              #corners <- data.frame("x" = c(1.13413, -3.56123, 157.593, 111.246, 1.13413), "y" = c(-0.0310049, 76.4206, -16.573, -63.6488, -0.0310049))
              corners <- data.frame("x" = c(borderLines$x,borderLines$x[1]) , "y" = c(borderLines$y,borderLines$y[1]))
              points(corners$x, corners$y, pch = 16, cex = 5, col = "red")
              lines(corners$x, corners$y, lty = 2, lwd = 2, col = "red")
            }
            
            
            treeListFile <- paste0(dirPath,basename(laz.path),"_trees_dbh.txt")
            if(file.exists(treeListFile)){
              trees <- read.table(treeListFile, header = TRUE)
              #trees$dCol <- trees$dbh/max(trees$dbh)
              trees$dCol <- rgb(50,250,50, maxColorValue = 255)
              trees$dCol[trees$id >= 9000] <- rgb(255,200,0, maxColorValue = 255)
              
              dbh.scaler <- 5
              if (greenTreeLocations) {
                points(trees$x, trees$y, pch = 16, cex = log(trees$dbh)/log(dbh.scaler)/1.5, col = trees$dCol)
              }
            }
            
            # central cross
            clip(x1 = -1, x2 = 1,
                 y1 = -2, y2 = 3)
            abline(v=0, h=0, lwd = 2)
            
          }, error = function(error_condition) 
          {
            cat("Error in creating this image, next loop!\n")
            next()
          })
        cat("\n")
        dev.off()
      }
    }
    cat("done")
    
    ### renaming and converting background image ####
    
    setStringApp <- paste0("_ct00_", -topLeftX*pixelUnit_cm, 
                           "_", topLeftY*pixelUnit_cm, 
                           "_", pixelUnit_cm*100, 
                           "_", modeApp)
    
    

    # copy first slice to main directory
    fromFile <- paste0(picPath.app.temp, outName, slices[1], ".png")
    
    # file.copy(fromFile, toFile, overwrite = T)
    
    if(createJPG){
      if(round(spanX_cm_drawn*pixelUnit_cm * spanY_cm_drawn*pixelUnit_cm / 1000000) > 600){
        # if image larger than 600 Mio pixel R fails to convert to jpg, better leave png manually!
        # this failure breaks down R session, so better be careful!
        cat(" - exceeding image limit, only renaming .png to .jpg (convert it yourself!)\n")
        toFile <- paste0(appPath,"/bgr_", tolower(outName), setStringApp, ".jpg")
        file.copy(fromFile, toFile, overwrite = T)
      } else {
        # read png and convert to jpg
        
        cat(" - converting .png to .jpg of quality", round(jpgQuality*100,1), "%...\n")
        toFile <- paste0(appPath,"bgr_", tolower(outName), setStringApp, ".jpg")
        img <- readPNG(fromFile)
        writeJPEG(img, target = toFile, quality = jpgQuality)
        rm(img)
        
      }
      
      
      gc()
    } else {
      if(createTIFF){
        cat(" - renaming .png to .tif file...\n")
        toFile <- paste0(appPath,"/bgr_", tolower(outName), setStringApp, ".tif")
      } else {
        cat(" - renaming .png file...\n")
        toFile <- paste0(appPath,"/bgr_", tolower(outName), setStringApp, ".png")
      }
      file.copy(fromFile, toFile, overwrite = T)
    }
    cat("Done in a ")
    print.difftime(round(Sys.time() - bgrTime, 1))
    cat("\n")
    unlink(fromFile)
    gc()
  }
  
  
  # WRITE COLORED LAZ SECTION ####
  if(writeColoredLAZ){
    
    
    taSet <- Sys.time()
    # synchronize these settings with the app
    tiles.spacingX <- 10
    tiles.spacingY <- 5
    tiles.points_per_meter <- 5
    
    
    
    
    coloredLasPath <- paste0(dirPath, "app/coloredLAS/")
    if(!exists(coloredLasPath)) dir.create(coloredLasPath)
    slices <- c("_clusterSlice_100to300.laz", "_clusterSlice_120to140.laz", "_clusterSlice_300to500.laz", "_ground.laz")
    {
      cat("\nColoring slices for manual tree finding:")
      cat(":\n")
    }
    
    minX <- 10000
    maxX <- -10000
    maxY <- -10000
    minY <- 10000
    
    for(i in 1:length(slices)){
      
      if(writeLASuncompressed){
        # uncompressed las file eats a lot of storage
        sliceOutFile <- paste0(fileFinder, substr(slices[i], 1, nchar(slices[i])-1), "s")
      } else {
        # compressed laz file for sending to another pc
        sliceOutFile <- paste0(fileFinder, substr(slices[i], 1, nchar(slices[i])-1), "z")
      }
      
      
      cat("    .* ", sliceOutFile, "... ", sep = "")
      
      
      
      
      ta <- Sys.time()
      laz.slice <- paste0(groundPath, fileFinder, slices[i])
      co <- capture.output(sliceLAS <- readLAS(laz.slice, select = selector))
      if(i == 1){
        allMaxInt <- max(sliceLAS@data$Intensity)
        minX <- sliceLAS@header@PHB$`Min X`
        maxX <- sliceLAS@header@PHB$`Max X`
        maxY <- sliceLAS@header@PHB$`Max Y`
        minY <- sliceLAS@header@PHB$`Min Y`
      } else {
        minX <- min(minX, sliceLAS@header@PHB$`Min X`)
        maxX <- max(maxX, sliceLAS@header@PHB$`Max X`)
        maxY <- max(maxY, sliceLAS@header@PHB$`Max Y`)
        minY <- min(minY, sliceLAS@header@PHB$`Min Y`)
      }
      
      if(slices[i] == "_clusterSlice_120to140.laz"){
        # color RED
        maxInt <- allMaxInt * 0.3
        sliceLAS@data$R <- 255L
        sliceLAS@data$G <- as.integer(255 - (sliceLAS@data$Intensity/maxInt) * 255)
        sliceLAS@data$G[sliceLAS@data$G < 0] <- 0L
        sliceLAS@data$B <- sliceLAS@data$G
      } else {
        # color intensity BLACK
        if(slices[i] == "_clusterSlice_100to300.laz"){
          maxInt <- allMaxInt * 0.5
          # regular intensity
        } else {
          # upper points (> 3 m) have less intensity, so make them darker (else look too whiteish)
          maxInt <- allMaxInt * 0.25
        }
        sliceLAS@data$G <- as.integer(255 - (sliceLAS@data$Intensity/maxInt) * 255)
        sliceLAS@data$G[sliceLAS@data$G < 0] <- 0L
        sliceLAS@data$R <- sliceLAS@data$G
        sliceLAS@data$B <- sliceLAS@data$G
        #sliceLAS@data$B <- 120L # blueish and yellow if enabled
      }
      sliceLAS@header@PHB$`Point Data Format ID` <- 3
      
      
      writeLAS(sliceLAS, paste0(coloredLasPath, sliceOutFile))
      te <- Sys.time()
      cat("done in", round(as.numeric(te-ta, units="secs"), 1), "sec!\n")
    }
    
    
    ### Tiles lines las
    {
      ta <- Sys.time()
      gridOutFile <- paste0(fileFinder, "_grid.las")
      outGrid.path <- paste0(coloredLasPath,gridOutFile)
      
      cat("    .* ", gridOutFile, "... ", sep = "")
      
      
      minX_l <- floor(minX/tiles.spacingX)*tiles.spacingX
      maxX_l <- ceiling(maxX/tiles.spacingX)*tiles.spacingX
      xRow <- seq(minX_l, maxX_l, by = tiles.spacingX)
      
      minY_l <- floor(minY/tiles.spacingY)*tiles.spacingY
      maxY_l <- ceiling(maxY/tiles.spacingY)*tiles.spacingY
      yRow <- seq(minY_l, maxY_l, by = tiles.spacingY)
      
      nTiles <- (length(xRow)-1)*(length(yRow)-1)
      cat(paste0("(", nTiles, " tiles) "))
      
      
      lines <- data.frame("x" = 1, 
                          "y" = rep(seq(minY_l, maxY_l, by = 1/tiles.points_per_meter), length(xRow)))
      lines$x <- (sort(rep(xRow, length(lines$y)/length(xRow))))
      
      lines2 <- data.frame("x" = rep(seq(minX_l, maxX_l, by = 1/tiles.points_per_meter), length(yRow)), "y" = 0)
      lines2$y <- (sort(rep(yRow, length(lines2$x)/length(yRow))))
      
      lines3 <- rbind(lines, lines2)
      #plot(lines3)
      colnames(lines3) <- c("X", "Y")
      
      
      
      # attaching z from DTM
      tryCatch(
        {
          
          dtm_z <- raster(paste0(groundPath, fileFinder, "_ground_min.grd"))
          lines3$Z <- round(extract(dtm_z, SpatialPoints(data.frame(x = lines3$X, y = lines3$Y))) + 1.3, 3)
          lines3$Z[is.na(lines3$Z)] <- 0
          }, error = function(error_condition) {
          cat("Error in reading the dtm-model, next loop!")
          return()
        })
      
      lines3$Intensity <- 1L
      lines4 <- lines3
      lines4$Z <- lines4$Z + 0.7 # 2m
      lines4$Intensity <- 2L
      lines5 <- lines4
      lines5$Z <- lines5$Z + 2
      lines5$Intensity <- 3L
      lines6 <- lines4
      lines6$Z <- lines4$Z - 2 # ground line
      lines6$Intensity <- 0L
      lines3 <- rbind(lines6, lines3, lines4, lines5)
      lines3 <- lines3[!is.na(lines3$Z),]
      
      
      
      co <- capture.output(gridLAS <- LAS(data = lines3), type = "message")
      #plot(gridLAS)
      writeLAS(gridLAS, outGrid.path)
      
      te <- Sys.time()
      cat("done in", round(as.numeric(te-ta, units="secs"), 1), "sec!\n")
      
    }
    
    teSet <- Sys.time()
    cat("         required in total", round(as.numeric(teSet-taSet, units="secs"), 1), "seconds.\n")
  }
  
  if(!is.na(fileFinder[1])){
    
  } else {
    
  }
  
  
  # CREATE TREE VOLUME TXT LIST SECTION ####
  {
    
    
    if(trees.path == ""){
      trees <- data.frame()
      trees_h_file <- NULL
      for(j in 1:length(fileFinder)){
        trees.file <- paste0(dirPath,fileFinder[j],"_ALLGO_100to300_dbh/trees_dbh.txt")
        if(!file.exists(trees.file)){
          trees.file <- paste0(dirPath,fileFinder[j],"/trees_dbh.txt")
          if(!file.exists(trees.file)){
            warning("Directory doesn't contain a trees_dbh.txt file!")
            cat("\n\n")
            return()
          }
        }
        check.file <- paste0(dirname(trees.file),"/trees_measured_species.txt")
        if(file.exists(check.file)){
          trees.file <- check.file
          cat("Working with species file!\n")
        } else {
          check.file <- paste0(dirname(trees.file),"/trees_measured.txt")
          if(file.exists(check.file)){
            trees.file <- check.file
            cat("Working with segmented crown file!\n")
          } else {
            trees_h_file <- list.files(dirname(trees.file), pattern = "trees_height", full.names = T)
          }
        }
        if(length(trees_h_file) > 0){
          if(j == 1) cat("Working with file", basename(trees_h_file[1]),"\n")
          trees <- rbind(trees, read.table(trees_h_file[1],
                                           header = T, dec = ".", sep = "\t"))
        } else {
          trees <- rbind(trees, read.table(trees.file, header = T, dec = "."))
        }
      }
      
    } else {
      
      
      if(endsWith(trees.path, "csv")){
        
        
        
        cat("Working with csv... ")
        trees <- read.csv(trees.path, header = T)
      } else {
        tryCatch(
          {
            # read in raster file
            cat("Working with table... ")
            trees <- read.table(trees.path, header = T, dec = ".", sep = "\t")
          }, error = function(error_condition) {
            cat("no - reading csv...")
            trees <- read.csv(trees.path, header = T)
          })
      }
      
      cat("done!\n")
      
      
    }
    
    
    if(is.element("species", colnames(trees))){
      if(typeof(trees$species) == "integer"){
        trees$species_text <- treeSpecies(trees$species)
        cat("Adding species_text column from species...\n")
      } else {
        trees$species_text <- trees$species
        trees$species <- as.integer(treeSpeciesNumber(trees$species_text))
        cat("Converting species number column...\n")
      }
    } else {
      warning("No tree species information - considering all as SPRUCE (1) Green.")
      
      trees$species <- 1
      trees$species_text <- "Fi"
    }
    trees$species[trees$species == 0] <- 1 # unknown species is spruce
    trees$species[trees$species == -1] <- 1 # unknown species is spruce
    #cat("NO - only Spruces!\n")
    #trees$species <- 1
    #trees$species_text <- "Fi"
    if(length(unique(trees$species)) > 1){
      print(table(trees$species_text))
      print(table(trees$species))
    }
    
    if(is.element("dbh_cm", colnames(trees))){
      trees$dbh <- trees$dbh_cm
    } 
    if(is.element("h_m", colnames(trees))){
      trees$est.height <- trees$h_m
    } 
    
    if(is.element("DBH", colnames(trees))){
      # corrections for the meta zenodo plots
      trees$dbh <- trees$DBH
      
      if(is.element("height", colnames(trees))){
        trees$est.height <- trees$height
        naHeights <- is.na(trees$est.height)
        if(sum(naHeights) > 0){
          trees$est.height[naHeights] <- 25
          cat("Taking tree height and correcting", sum(naHeights), "trees to 25 m!\n")
        }
      }
      
      naTrees <- is.na(trees$treeID)
      trees$treeID[naTrees] <- c(1:sum(naTrees))+4000
      trees$dbh[naTrees] <- 4.4
      trees$est.height[naTrees] <- 5
      
      if(is.element("treeID", colnames(trees))){
        trees$id <- trees$treeID
      }
    }
    
    if(!is.element("z", colnames(trees))){
      cat("Attaching missing z-values from ground model.\n")
      if(!exists("dtm")) dtm <- raster(groundModel.path)
      trees$z <- round(extract(dtm, SpatialPoints(data.frame(x = trees$x, y = trees$y))) + 1.3, 3)
      
      # fix na of z values
      naZ <- which(is.na(trees$z))
      if(length(naZ > 0)){
        for(i in seq_along(naZ)){
          treePos <- data.frame(x = trees$x[naZ[i]], y = trees$y[naZ[i]])
          radius <- c(0.3, 0.5)  # 30 and 50 cm circle
          n_points <- 72  
          # Create angles evenly spaced around the circle
          angles <- seq(0, 2*pi, length.out = n_points + 1)[- (n_points + 1)]  # exclude last point to avoid duplication
          # Compute x and y coordinates for each point
          circle_points <- data.frame(
            x = treePos$x + radius * cos(angles),
            y = treePos$y + radius * sin(angles)
          )
          #plot(circle_points, asp = 1)
          #text(x = circle_points$x, y = circle_points$y, labels = newZ)
          outZ <- round(extract(dtm, SpatialPoints(circle_points)) + 1.3, 3)
          averagedZ <- mean(newZ, na.rm = T)
          if(is.na(averagedZ)){
            radius <- c(1)  # 1 m circle
            n_points <- 36  
            angles <- seq(0, 2*pi, length.out = n_points + 1)[- (n_points + 1)]  # exclude last point to avoid duplication
            circle_points <- data.frame(
              x = treePos$x + radius * cos(angles),
              y = treePos$y + radius * sin(angles)
            )
            #plot(circle_points, asp = 1)
            #text(x = circle_points$x, y = circle_points$y, labels = newZ)
            outZ <- round(extract(dtm, SpatialPoints(circle_points)) + 1.3, 3)
            averagedZ <- mean(newZ, na.rm = T)
            if(is.na(averagedZ)){
              averagedZ <- -100
            }
          }
          trees$z[naZ[i]] <- averagedZ
        }
      }
      }
    trees$z[is.na(trees$z)] <- 0
    
    if(is.element("randomCol", colnames(trees))){
      trees$comment <- trees$randomCol
    }
    if(!is.element("comment", colnames(trees))){
      warning("Adding empty field for comment")
      trees$comment <- ""
    }
    
    
    if(is.element("BHD", colnames(trees))){
      warning("Converting BHD to DBH - is input file in German?")
      trees$dbh <- trees$BHD
    }
    
    #safeT <<- trees
    if(!is.element("est.height", colnames(trees))){
      warning("No tree height information - considering all as 25 m high.")
      trees$est.height <- 25
    }
    
    if(!is.element("z", colnames(trees))){
      warning("No z-coordinate - considering all at 0 m Altitude.")
      trees$z <- 0
    }
    
    #safeT <<- trees
    # compute volume
    trees$fz <- fPoll(trees$species, trees$dbh, trees$est.height)
    #safeT <<- trees
    trees$v <- (trees$dbh/100)^2 * pi / 4 * trees$est.height * trees$fz
    #safeT <<- trees
    trees$v_denzin <- trees$dbh^2/1000
    #safeT <<- trees
    
    
    
    plot(trees$est.height ~ trees$dbh)
    #identify(trees$dbh, trees$est.height) #pick single outliers
    if(1 == 2){
      #install.packages("gatepoints")
      #library(gatepoints) # for selecting outliers in the dbh/h plot
      X11() #open new plot window
      selectedPoints <- fhs(trees) #draw a polygon to find outliers
      selectedPoints #outliers
    }
    
    
    
    # if(!is.element("inside", colnames(trees))){
    if(TRUE){
      if(circleRadius > 0){
        area_m <- 0
        cat("Checking inside trees for", circleRadius, "m radius...")
        trees$dist <- sqrt(trees$x^2 + trees$y^2)
        trees$inside <- F
        trees$inside[trees$dist <= circleRadius] <- T
        
        area_m <- round(circleRadius^2*pi,2)
        
      } else {
        
        if(!is.na(fileFinder[1])){
          if(fileFinder[1] == "b_Altmuenster" || fileFinder[1] == "a_Rosalia" || fileFinder[1] == "c_Kleinarl" || fileFinder[1] == "d_Brixental"){
            cat("Setting our borders for the trees")
            trees$inside <- inside.owin(trees$x, trees$y, trees_edges_buffer)
          } else {
            cat("No border lines found, setting all as inside!\n")
            trees$inside <- T
          }
        } else {
          trees$inside <- T
        }
      }
    }
    
    plot(trees$y ~ trees$x, asp = 1)
    trees$inside <- trees$inside == 1 #convert to boolean
    points(trees$x[trees$inside], trees$y[trees$inside], col = "red")
    trees_out <- trees[!trees$inside,]
    trees_in <- trees[trees$inside,]
    points(trees_out$y ~ trees_out$x, col = "blue", pch = 13)
    if(exists("area_m")){
      title(paste0(fileFinder, " ", area_m, "m2    trees:", length(trees_in$x), "in  ", length(trees_out$x), "out"))
    } else {
      title(paste0(fileFinder, "    trees:", length(trees_in$x), "in  ", length(trees_out$x), "out"))
      
    }
    
    
    
    if(!is.element("selected", colnames(trees))){
      cat("No prior information - setting all trees as unselected!\n")
      trees$selected <- FALSE
    }
    
    if(2==1){
      # reading a list and setting trees pre-selected
      sel <- read.table("F:/sel.txt", head = T)
      trees$selected[is.element(trees$id, sel$x)] <- TRUE
      table(trees$selected)
    }
    
    #plot(trees$v ~ trees$dbh)
    
    total_volume <- sum(trees_in$v)
    total_volume_denzin <- sum(trees_in$v_denzin)
    cat("Total volume in 1 ha is", total_volume, "Vfm - Denzin volume is", total_volume_denzin, "Vfm.\n")
    try(hist(trees_in$dbh))
    
    
    sum(trees_in$v)
    shortList <- data.frame("id" = trees_in$id,
                            "x" = round(trees_in$x,3),
                            "y" = round(trees_in$y,3),
                            "dbh" = trees_in$dbh,
                            "h" = trees_in$est.height,
                            "vol" = round(trees_in$v,3),
                            "fz" = round(trees_in$fz,3),
                            "hd" = round(trees_in$est.height*100/trees_in$dbh,2))
    
    # write.table(shortList, file = paste0(dirPath,"trees_",filefinder,"_inside.txt"),
    #             row.names = FALSE, sep = "\t")
    hist(shortList$hd)
    sum(shortList$v) #759.211Vfm
    table(trees_in$species_text) #601Bu 591Fi, 34La, 13Ta
    length(shortList$v) #1239 inside
    mean(shortList$hd) #111.24
    
    appList <- data.frame("x" = trees$x,
                          "y" = trees$y,
                          "z" = trees$z,
                          "id" = trees$id,
                          "comment" = trees$comment,
                          "species" = trees$species_text,
                          "dbh" = trees$dbh,
                          "height" = trees$est.height,
                          "vol" = round(trees$v,3),
                          "fz" = round(trees$fz,3),
                          "inside" = as.integer(trees$inside),
                          "selected" = as.integer(trees$selected))
    
    if(eraseSpecies) {
      appList$species <- "NA"
      appList$comment <- ""
    }
    
    if(highlightNew){
      appList$species[appList$id >= 9000] <- "GX"
    }
    
    if(setSpecies != ""){
      appList$species <- setSpecies
    }
    
    
    
    write.table(appList, file = paste0(dirPath,"app/trees_",tolower(outName),".txt"),
                row.names = FALSE, sep = "\t")
    cat("  -> set", outName, "done!")
   
    if(length(list.files(picPath.app.temp)) == 0) {
      unlink(picPath.app.temp, recursive = T)
    } 
    if(wait5s){
      Sys.sleep(2)
      cat("\n")
      Sys.sleep(2)
      cat("\n")
      Sys.sleep(1)
    }
    ############################################
  }
  cat("\n")
}




#
# opens multiple height slices in Cloudcompare
# for better visualisation and finding omitted trees
#
# opens following slices:
#   * red DBH (120 - 140 cm)
#   * bw 100 - 300 cm
#   * bw 300 - 500 cm
#  (* bw 0 - 100 cm (ground))
#   ' grid 5 x 5 m (5 x 10 m)
# 
# optionally opens ground (0 - 100 cm) by setting    withGround <- T
#' @export
openColoredLAZ <- function(fileFinder, 
                           withGround = T, 
                           dirPath = paste0(getwd(), "/")){
  #library(lidR)
  library(PBSmodelling) #open files
  
  # # # open and explore # # #
  # the path for ground and vegetation files
  
  #appPath <- paste0(inPath, "trees/")
  appPath <- paste0(dirPath, "app/")
  coloredLasPath <- paste0(dirPath, "app/coloredLAS/")
  
  
  
  # DEFINE HERE WHICH FILE TO OPEN
  #EB openThis <- 159
  #EB fileFinder <- paste0("eb", sprintf("%03d", openThis))
  #EB fileFinder <- fileFinders[openThis]
  #EB allFileFinderFiles <- list.files(path = appPath, pattern = fileFinder)
  allFileFinderFiles <- list.files(path = appPath, pattern = paste0(tolower(fileFinder), ".txt"))
  if(length(allFileFinderFiles) > 0){
    if(length(allFileFinderFiles) > 1){
      allFileFinderFiles <- allFileFinderFiles[order(allFileFinderFiles, decreasing = T)]
      treeListPath <- paste0(appPath, allFileFinderFiles[1])
    } else {
      treeListPath <- paste0(appPath, allFileFinderFiles)
    }
    cat("Loading tree list file:", treeListPath)
    if(withGround){
      cat(" with ground.\n")
    } else {
      cat(" without ground.\n")
    }
  } else {
    treeListPath <- ""
  }
  
  
  las.slice.cluster.1 <- paste0(coloredLasPath, fileFinder, "_clusterSlice_120to140.las")
  las.slice.cluster.2 <- paste0(coloredLasPath, fileFinder, "_clusterSlice_100to300.las")
  las.slice.cluster.3 <- paste0(coloredLasPath, fileFinder, "_clusterSlice_300to500.las")
  las.ground <- paste0(coloredLasPath, fileFinder, "_ground.las")
  if(!file.exists(las.slice.cluster.2)){
    las.slice.cluster.1 <- paste0(coloredLasPath, fileFinder, "_clusterSlice_120to140.laz")
    las.slice.cluster.2 <- paste0(coloredLasPath, fileFinder, "_clusterSlice_100to300.laz")
    las.slice.cluster.3 <- paste0(coloredLasPath, fileFinder, "_clusterSlice_300to500.laz")
    las.ground <- paste0(coloredLasPath, fileFinder, "_ground.laz")
  }
  txt.trees.dbh <- treeListPath
  # openFile(las.slice.cluster.slope, las.slice.raw)
  
  
  #open cloudcompare slice file with tree list
  #openFile(paste0(dirPath, fileFinder, "_ALLGO_100to300_dbh/"))
  
  if(withGround){
    lasOpen <- (c(las.slice.cluster.1, las.slice.cluster.2, las.slice.cluster.3, 
                  las.ground, txt.trees.dbh))
  } else {
    lasOpen <- (c(las.slice.cluster.1, las.slice.cluster.2, las.slice.cluster.3, 
                  txt.trees.dbh))
  }
  
  file.exists(lasOpen)
  
  gridFile.path <- paste0(coloredLasPath, fileFinder, "_grid.las")
  if(file.exists(gridFile.path)){
    lasOpen <- c(lasOpen, gridFile.path)
  }
  
  
  File.exe <- "Cloudcompare.exe"
  Path.CloudCompare <- "\"C:/Program Files/CloudCompare/\""
  Path.exe <- paste0(Path.CloudCompare, File.exe)
  (command.exe <- paste0(" ", Path.exe, " ", 
                         gsub(",","",toString(lasOpen))))
  
  cat(command.exe)
  system(command.exe, wait = FALSE)
  
  
}




#' 
#' Creates a list of how many new trees were added ( id > 9000 )
#' and also gives the mean, max and median of new and old DBH
#' 
#' @param dir_completedInputLists path to the list of exported trees (new trees > 9000)
#' @export
statsCompletedTrees <- function(dir_completedInputLists){
  
  library(dplyr)
  timeAn1 <- Sys.time()
  cat("Analyzing statistics of completed trees in directory ",dir_completedInputLists,"\n")
  
  treeLists <- list.files(dir_completedInputLists, pattern = "_db.txt", full.names = T)

  fileFinders <- sub("_.*", "", basename(treeLists))
  
  date_time <- sub(".*_.*_.*_([0-9]{8}_[0-9]{6}).*", "\\1", 
                   basename(treeLists))
  
  inputFiles <- data.frame(
    fileFinders_list = fileFinders,
    date_time = date_time,
    listPath = treeLists,
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
  cat(" (of total", length(treeLists), "lists)\n")
  
  
  fileFinders <- uniqueFiles_total$fileFinders_list
  treeLists <- uniqueFiles_total$listPath
  
  completedTrees <- data.frame("fileFinder" = fileFinders, 
                               "newTrees" = 0L, 
                               "meanDBH" = 0, "maxDBH" = 0, "medianDBH" = 0, 
                               "oldTrees" = 0L, 
                               "meanDBHOlds" = 0, "maxDBHOlds" = 0, "medianDBHOlds" = 0)
  
  
  for(i in seq_along(treeLists)){
    nowTrees <- read.table(treeLists[i], header = T, sep = "\t")
    
    newTrees <- nowTrees[nowTrees$id >= 9000,]
    oldTrees <- nowTrees[nowTrees$id < 9000,]
    
    
    if(nrow(newTrees) > 0){
      completedTrees$newTrees[i] <- nrow(newTrees)
      completedTrees$meanDBH[i] <- mean(newTrees$dbh)
      completedTrees$medianDBH[i] <- median(newTrees$dbh)
      completedTrees$maxDBH[i] <- max(newTrees$dbh)
    } else {
      cat(paste0("No new trees in set ", fileFinders[i], "!\n", sep = ""))
    }
    
    if(nrow(oldTrees) > 0){
      completedTrees$oldTrees[i] <- nrow(oldTrees)
      completedTrees$meanDBHOlds[i] <- mean(oldTrees$dbh)
      completedTrees$medianDBHOlds[i] <- median(oldTrees$dbh)
      completedTrees$maxDBHOlds[i] <- max(oldTrees$dbh)
    } else {
      cat(paste0("No OLD trees automatically found in ", fileFinders[i], "!\n", sep = ""))
    }
  }
  
  meanLine <- colMeans(completedTrees[, c(2:ncol(completedTrees))])
  sumLine <- colSums(completedTrees[, c(2:ncol(completedTrees))])
  headerLines <- rbind(c(123, sumLine),
                       c(123, meanLine))
  headerLines <- data.frame(headerLines)
  colnames(headerLines)[1] <- "fileFinder"
  headerLines$fileFinder <- c("sum", "mean")
  headerLines[1, c(3:5, 7:9)] <- NA
  dout <-rbind(headerLines, completedTrees)
  
  cat("\n\nAll sets analyzed.\n")
  cat("In total there were", sum(completedTrees$newTrees), "new trees added\n")
  cat(paste0("   (", round(sum(completedTrees$newTrees) / nrow(completedTrees), 1), " per plot)\n\n"))
  
  write.csv2(dout, paste0(dir_completedInputLists, 
                          "../completed_tree_statistics_", 
                          format(Sys.time(), "%Y%m%d_%H%M"), ".csv"), 
             row.names = F
  )
  
  timeAn2 <- Sys.time()
  cat("Analysis done in a ")
  print.difftime(round(timeAn2 - timeAn1), 1)
  cat("\n\n")
  
}



#cat("Succesfully updated function createAppFiles() V25.01 with clip window!\n")
