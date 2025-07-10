

if(!exists("sliVox")){
  sliVox <<- NA # holds the voxelized slice, to extract per cluster points for dbh measuring faster
}


#' @export
removeUmlaut <- function(inputString){
  out <- stringi::stri_replace_all_fixed(
    inputString,
    c("\U00E4", "\U00F6", "\U00FC", "\U00C4", "\U00D6", "\U00DC", "\U00DF"),
    c("ae", "oe", "ue", "Ae", "Oe", "Ue", "ss"),
    vectorize_all = FALSE)
  if(out != inputString){
    warning("Modified the fileFinder. There should be no special characters!\n")
  }
  return(out)
}

#' Individual tree detection
#'
#' Reads ground and vegetation point cloud and performs
#' tree localisation by density cluster analysis and
#' produces trees_dbh.txt file
#'
#' @param allDBHs add all measured diameters (gam, circle, ellipsis) to the output file trees_dbh.txt
#* @param nr_cores how many cores to use for parallel cluster analysis
#' @export
#' @useDynLib edc
clustSplit <- function(fileFinder, allDBHs = FALSE, allFiles = FALSE,
                       clipHeight = 3, bottomCut = 1, ipad = FALSE, nr_cores = 0,
                       bushPreparation = FALSE, filterSOR = FALSE, filterINT = 0, ref = NA, ref.plot_id = NA,
                       cutWindow = c(-1000, -1000, 2000), numberOfPoints = 300, heightExtent = 1.3, TLS = FALSE,
                       silent = TRUE, fast = TRUE, retainPointClouds = TRUE, dirPath = paste0(getwd(), "/")){



  
  # only for debugging
  if(FALSE){
    allDBHs = FALSE
    allFiles = FALSE
    
    clipHeight = 3
    bottomCut = 1
    ipad = FALSE
    nr_cores = 0
    
    bushPreparation = FALSE
    filterSOR = FALSE
    filterINT = 0
    ref = NA
    ref.plot_id = NA
    
    cutWindow = c(-1000, -1000, 2000)
    numberOfPoints = 300
    heightExtent = 1.3
    TLS = FALSE
    
    silent = TRUE
    fast = TRUE
    retainPointClouds = TRUE
    dirPath = paste0(getwd(), "/")
  }

  if(ipad){
    clipHeight <- 2.2
    bottomCut <- 0.2
    heightExtent <- 1.0
    numberOfPoints <- 0
    filterINT <- 0
  }

  gloStart <- Sys.time()
  fileFinder <- removeUmlaut(fileFinder)

  #cat("Refering job for total", length(trees.file), "cases...\n")
  referenced <- FALSE
  if(!is.na(ref)) referenced <- TRUE
  setStr <- generateSetString(fileFinder = fileFinder, mode = "ALLGO",
                              clipHeight = clipHeight, bottomCut = bottomCut,
                              bushPreparation = bushPreparation,
                              filterSOR = filterSOR, cutWindow = cutWindow, silent = TRUE)
  dbhPath <- paste0(dirPath, setStr, "_dbh/")
  if(!dir.exists(dbhPath)) dir.create(dbhPath)
  sink(paste0(dbhPath, fileFinder, "_clustSplit", format(Sys.time(), "%Y%m%d_%H%M"), "_Rcons.txt"), append = TRUE, split = TRUE)

  cat("\nStarting ALLGO Stem Detection and DBH Estimation for", fileFinder, "\n")
  cat("Today is", format(Sys.time()), "\n")
  #cat("\nLet's go with ALLGO cluster splitting :D!\n")
  cat("Working path is ", dbhPath, "\n")
  if(!exists("LAS_veg")){
    LAS_veg <<- NA
    LAS_veg_name <<- "blank"
  }



  {
    #cat("Creating folder structure... ")

    path.output.cluster.end <- paste0(dbhPath, "fineCluster/")
    path.output.cluster.endgraph <- paste0(dbhPath, "graphSlice/")

    if(!dir.exists(path.output.cluster.end)) dir.create(path.output.cluster.end)
    if(allFiles){
      if(!dir.exists(paste0(dbhPath, "residuals/"))) dir.create(paste0(dbhPath, "residuals/"))
      if(!dir.exists(path.output.cluster.endgraph)) dir.create(path.output.cluster.endgraph)
    }

    if(referenced){
      plotReferenceFile(fileFinder = fileFinder, ref.file = ref, ref.plot_id = ref.plot_id,
                        cutWindow = cutWindow, writePNG = TRUE, pathPNG = dbhPath)
    }
    #cat("done!\n")
  }


  if(file.exists(paste0(dbhPath, "slice_cluster.laz"))){
    cat("Skipping roughCluster(), loading old slice_cluster.laz file... ")
    sliVox <<- readLAS(paste0(dbhPath, "slice_cluster.laz"))
    cat("done!\n")
  } else {
    roughCluster(fileFinder, dbhPath = dbhPath, ipad = ipad, allFiles = allFiles,
                 clipHeight = clipHeight, bottomCut = bottomCut, referenced = referenced,
                 numberOfPoints = numberOfPoints, heightExtent = heightExtent, TLS = TLS,
                 bushPreparation = bushPreparation, filterSOR = filterSOR, filterINT = filterINT,
                 cutWindow = cutWindow, silent = silent, retainPointClouds = retainPointClouds, dirPath = dirPath)
  }

  diameterBeast(fileFinder, dbhPath = dbhPath, referenced = referenced, ipad = ipad, allFiles = allFiles,
                bushPreparation = bushPreparation, filterSOR = filterSOR, nr_cores = nr_cores,
                cutWindow = cutWindow, silent = silent, fast = fast, dirPath = dirPath)

  fineCluster(fileFinder, dbhPath = dbhPath, allDBHs = allDBHs, referenced = referenced, 
              bushPreparation = bushPreparation, filterSOR = filterSOR, nr_cores = nr_cores, 
              cutWindow = cutWindow, silent = silent, dirPath = dirPath)

  if(!retainPointClouds){
    LAS_veg <<- NA #removing old LAS_veg to ensure performance
    gc()
  }

  gloStop <- Sys.time()
  cat("\nTree detection completed in a ")
  print.difftime(round(gloStop - gloStart, 1))

  if(referenced){
    cat("Referencing not implemented yet! \n\n")
  }

  cat("\n\n\n")
  sink()
}


####################################################################################################
# 1 FIRST CLUSTERING ###############################################################################
####################################################################################################

roughCluster <- function(fileFinder, dbhPath, ipad = FALSE, allFiles = FALSE,
                         clipHeight = 3, bottomCut = 1, numberOfPoints = 300, heightExtent = 1.3,
                         referenced = FALSE, bushPreparation = FALSE, TLS = FALSE, filterSOR = FALSE, filterINT = 0,
                         fast = TRUE,
                         cutWindow = c(-1000, -1000, 2000), silent = FALSE, retainPointClouds = FALSE, dirPath = paste0(getwd(), "/")){

  start <- Sys.time()
  cat("\nStarting roughCluster()\n")

  setStr <- generateSetString(fileFinder = fileFinder, mode = "ALLGO",
                              clipHeight = clipHeight, bottomCut = bottomCut,
                              bushPreparation = bushPreparation,
                              filterSOR = filterSOR, cutWindow = cutWindow, silent = TRUE)



  ### READING or CREATING SLICE FILE FROM VEGETATION ####

  start <- Sys.time()
  groundPath <- v.env$groundPath
  slicePath <- paste0(dirPath, groundPath, fileFinder, "_clusterSlice_", bottomCut*100, "to", clipHeight*100, ".laz")
  slicePathAlt <- paste0(dirPath, groundPath, fileFinder, "_clusterSlice_", bottomCut*100, "to", clipHeight*100, "ORI.laz")
  if(file.exists(slicePathAlt)) slicePath <- slicePathAlt
  voxSlicePath <- paste0(dirPath, groundPath, fileFinder, "_clusterSlice_", bottomCut*100, "to", clipHeight*100, "_vox.laz")
  voxSlicePath.slope <- paste0(dirPath, groundPath, fileFinder, "_clusterSlice_", bottomCut*100, "to", clipHeight*100, "_vox_slope.laz")

  cat("Checking for old slice file at", slicePath, "...\n")


  if(file.exists(slicePath)){
    cat("-> We use prevailing raw slice file... ")
    if(sum(cutWindow == c(-1000, -1000, 2000)) == 3){
      #file.copy(slicePath, paste0(dbhPath, "slice_raw.laz"))
      slice <- readLAS(slicePath, select = "xyzcit0")
      cat("complete file copied.\n")
    } else {
      cat("reading for clipping... ")
      slice <- readLAS(slicePath, select = "xyzcit0")
      numP <- slice@header@PHB$`Number of point records`
      cat("found originally", thMk(numP), "points.\n")
      cat("Clipping to dimensions", cutWindow, "leaves ")

      XL <- cutWindow[1]
      YL <- cutWindow[2]
      width <- cutWindow[3]
      slice <- filter_poi(slice, X > XL, X < XL + width, Y > YL, Y < YL + width)

      numP <- slice@header@PHB$`Number of point records`
      cat(thMk(numP), "points, writing out slice_raw.laz... ")
      writeLAS(slice, paste0(dbhPath, "slice_raw.laz"))
      cat("done!\n")
      gc()
    }
  }



  # if(file.exists(voxSlicePath) && !retainPointClouds){ #CHANGE HERE to always load file when RETAINPOINTCLOUDS = FALSE
  if(file.exists(voxSlicePath) && filterINT == 0){ #old one, please change! AT 21-05-18
    cat("-> We also use processed voxel cluster file, reading in... ")
    sliVox <<- readLAS(voxSlicePath, select = "xyzcit0")
    numP <- sliVox@header@PHB$`Number of point records`
    cat("done.\n")
    try(file.copy(voxSlicePath.slope, paste0(dbhPath, "slice_cluster_slope.laz")))
    cat("Found", length(unique(sliVox@data$cluster)), "clusters in", thMk(numP),
        "points, z-values ranging from", round(sliVox@header@PHB$`Min Z`, 2), "to", round(sliVox@header@PHB$`Max Z`, 2), "m.\n")
  } else {

    # Creating the basic slice file
    if(!file.exists(slicePath)){
      cat("-> We create new slice file! \nBy ")
      if(is.na(LAS_veg) || fileFinder != LAS_veg_name){
        cat("reading in: ", paste0(fileFinder, "_raw_veg.laz"), "...\n", sep = "")
        vegetation <- readLAS(file = paste0(dirPath, groundPath, fileFinder, "_raw_veg.laz"), select = "xyzcit0")

        if(retainPointClouds){
          cat("Retaining LAS_veg variable for ")
          LAS_veg <<- vegetation #retain big LAS file in memory, less performant
          LAS_veg_name <<- fileFinder
          cat(LAS_veg_name, "\n")
        }

        if(bottomCut < 1){
          cat("Also reading in: ", paste0(fileFinder, "_ground.laz"), "...\n", sep = "")
          ground <- readLAS(file = paste0(dirPath, groundPath, fileFinder, "_ground.laz"), select = "xyzcit0")
          vegetation <- rbind(vegetation, ground)
          rm(ground)
          gc()
        }



      } else {
        cat("using the prevailing \"LAS_veg\" variable to extract the", LAS_veg_name, "point cloud...\n")
        vegetation <- LAS_veg
      }



      tryCatch(
        {
          # read in raster file
          dtm_z <- raster(paste0(dirPath, groundPath, fileFinder, "_ground_min.grd"))
        }, error = function(error_condition) {
          cat("Error in reading the dtm-model, file not found!")
          return()
        })

      useCoarseGrid <- FALSE
      # NORMALIZATION
      tryCatch(
        {
          cat("Normalizing height... ")
          vegetation <- normalize_height(vegetation, dtm_z, na.rm = TRUE)
          cat("done!\n")
        }, error = function(error_condition) {
          useCoarseGrid <<- TRUE
          cat("-> problem with the normalisation, the minimal DTM model seems to be faulty...\n")
        })
      if(useCoarseGrid){
        # read in raster file 2, because first one produces some NAs at normalization
        cat("\"Roughly\" normalizing height - UNEXCEPTIONAL CASE!!! - (using the coarse 1 x 1 m dtm)...")
        tryCatch(
          {
            # reading raster file
            dtm_a <- raster(paste0(dirPath, groundPath, fileFinder, "_ground_rough.grd"))
          }, error = function(error_condition) {
            cat("Error in reading the rough dtm-model, file not found!")
            return()
          })

        tryCatch(
          {
            vegetation <- normalize_height(vegetation, dtm_a, na.rm = TRUE) # need to save it in that intermediate object or it cannot unnormalize anymore
            cat("done!\n")
          }, error = function(error_condition) {
            cat("Error in creating the rough dtm model!")
            return()
          })
      }

      #rm(dtm, dtm_c, dtm3)
      slice <- filter_poi(vegetation, Z < clipHeight, Z > bottomCut)
      slice_un <- unnormalize_height(slice)
      cat("Creating global slice file at", slicePath, "... ")
      writeLAS(slice_un, slicePath)
      cat("done.\n")

      numP <- slice@header@PHB$`Number of point records`

      if(sum(cutWindow == c(-1000, -1000, 2000)) == 3){
        cat("No cutting because of no cutWindow setting!\n")
      } else {
        cat("Just for output: clipping raw", thMk(numP), "points to dimensions", cutWindow, "leaves")
        XL <- cutWindow[1]
        YL <- cutWindow[2]
        width <- cutWindow[3]
        slice_un <- filter_poi(slice_un, X > XL, X < XL + width, Y > YL, Y < YL + width)
        numP <- slice_un@header@PHB$`Number of point records`
        cat(thMk(numP), "points.\n")
      }

      cat("Creating local cropped slice file at", slicePath, "... ")
      writeLAS(slice_un, paste0(dbhPath, "slice_raw.laz"))
      cat("done, now we proceed with full cloud again.\n")

      slice <- slice_un
      rm(vegetation)
      gc()

      stop <- Sys.time()
      print.difftime(round(stop - start, 1))
    }


    if(filterSOR){
      cat("Applying noise filter from inside point cloud (no separate settings specified):")
      slice <- filter_poi(slice, Classification < 2) #0 = never classified, #1 = unclassified (vegetation) #18 = noise
      cat("remaining", thMk(slice@header@PHB$`Number of point records`), "points (approx.", round(slice@header@PHB$`Number of point records`/numP*100, 1), "%).\n")
    }

    if(filterINT != 0){
      cat("Applying intensity filter for", filterINT, "percentile: \n")
      threshold <- quantile(slice@data$Intensity, 1 - filterINT/100) # all above are 5 %
      cat("Keeping", round(filterINT, 1), "% of all points with an intensity higher than", thMk(threshold), "...\n")
      nPointsOld <- slice@header@PHB$`Number of point records`
      slice <- filter_poi(slice, Intensity > threshold)
      nPointsNew <- slice@header@PHB$`Number of point records`
      percentRem <- nPointsNew/nPointsOld
      cat("There are", thMk(nPointsNew), "points remaining (equals only",
          round(percentRem*100, 1), "% of original data).\n")
    }

    ### VOXEL SYSTEMATICS ####

    vox.size <- 0.015
    if(ipad) {
      vox.size <- 0.01
    }

    cat("\nVoxelize points with a raster of", vox.size*100, "cm... ")
    {
      #t1 <- Sys.time()
      #thin4 <- voxelize_points(slice, vox.size)
      t2 <- Sys.time()
      #cat("lidR done")
      #print.difftime(t2-t1)
      thin5 <- tlsSample(slice, smp.voxelize(vox.size)) #bei TLS 0.015 # bei PLS ist 0.02 anders # auch da ist das mit smp.voxelize neu # iPAd 0.015
      t3 <- Sys.time()
      cat("done by treeLS. ")
      print.difftime(round(t3-t2,1))
    }
    sliVox <<- thin5
    rm(slice, thin5)
    gc()




    ### CLUSTERING STEP ####
    #Ordering points to identify the clustering structure
    cat("Identifying clusters with OPTICS alogrithm... ")
    t1 <- Sys.time()
    if(ipad){
      res <- optics(cbind(sliVox@data$X, sliVox@data$Y), eps = 0.03,  minPts = 90) #eps = 0.025,  minPts = 90
    } else {
      res <- optics(cbind(sliVox@data$X, sliVox@data$Y), eps = 0.025,  minPts = 90) #eps = 0.025,  minPts = 90
    }
    t2 <- Sys.time()
    cat("done. ")
    print.difftime(round(t2-t1, 1))
    print(res)

    #Cluster herausgeben
    if(ipad){
      res <- extractDBSCAN(res, eps_cl = .030)
    } else {
      res <- extractDBSCAN(res, eps_cl = .025) #0.025
    }




    cat("\nIn total there were ")
    # add information to which cluster point belongs
    sliVox <<- add_lasattribute(sliVox, x = res$cluster, name = "cluster", desc = "ID of first stem clusters")
    numClustBefore <- length(unique(sliVox@data$cluster))
    cat(numClustBefore, "clusters detected, they are filtered now to number of Points and z-extent:\n")
    #lidR::plot(sliVox, color = "cluster")


    # filtering height and numberpoints
    {
      cat("Reducing", thMk(sliVox@header@PHB$`Number of point records`), "points in slice")
      hilf.tab <- as.data.frame(table(sliVox@data$cluster), stringsAsFactors = FALSE)
      colnames(hilf.tab) <- c("cluster", "anz")
      hilf.tab$zRange <- aggregate(sliVox@data$Z, by=list(sliVox@data$cluster), FUN=function(x) max(x)-min(x))[, 2]
      hilf.tab <- hilf.tab[hilf.tab$cluster!=0, ]
      try(hilf.tab$cluster <- strtoi(hilf.tab$cluster))
      hilf.tab$keepLarger <- hilf.tab$anz>=numberOfPoints  # filtering clusters that have less than numberOfPoints points
      hilf.tab$keepZExtent <- hilf.tab$zRange >= heightExtent # filtering clusters that are spanning less than 1.30 in height
      hilf.tab$keep <- hilf.tab$keepLarger & hilf.tab$keepZExtent
      hilf.tab

      sliVox <<- filter_poi(sliVox, cluster%in%hilf.tab[hilf.tab$keep, ]$cluster)
      numClustAfter <- length(unique(sliVox@data$cluster))
      cat(" to", thMk(sliVox@header@PHB$`Number of point records`), "points, lost", numClustBefore - numClustAfter, "clusters.\n")

      cat("Creating global output file *_vox.laz containing", numClustAfter, "clusters... ")
      writeLAS(sliVox, voxSlicePath)
      cat("done.\n")
    }


  }



  if(sum(cutWindow == c(-1000, -1000, 2000)) == 3){
    sliVox <- sliVox
    numClustCut <- length(unique(sliVox@data$cluster))
    cat("No cutting because of no cutWindow setting, all", numClustCut, "clusters remain!\n")
  } else {
    XL <- cutWindow[1]
    YL <- cutWindow[2]
    width <- cutWindow[3]
    cat("Cutting to dimensions", cutWindow, "... ")
    sliVox <- filter_poi(sliVox, X > XL, X < XL + width, Y > YL, Y < YL + width)

    numClustCut <- length(unique(sliVox@data$cluster))
    cat("done, ", numClustCut, "clusters remain.\n\n")
  }



  cat("Drawing image cluster_all.png... ")
  {
    png(paste0(dbhPath, "cluster_all.png"), height = 4000, width = 4000)
    plot(sliVox@data$X, sliVox@data$Y, col=sliVox@data$cluster, asp=1, cex=0.5, cex.axis=4) #Zeichnen der geclusterten Punkte
    #Dazuschreiben der Clusternummer
    hilf <- aggregate(cbind(sliVox@data$X, sliVox@data$Y), by=list(sliVox@data$cluster), FUN=mean)
    colnames(hilf) <- c("Cluster", "x", "y")
    points(hilf$x, hilf$y,  col=hilf$Cluster, pch=8, cex=0.8)
    text(hilf$x, hilf$y,  labels=hilf$Cluster, cex=1.5, pos=4, offset=0.5)
    dev.off()
    cat("done!\n")
  }

  # NORMALIZATION
  cat("Normalizing finally height of the cluster slice again... ")
  tryCatch(
    {
      # read in raster file
      dtm_z <- raster(paste0(dirPath, groundPath, fileFinder, "_ground_min.grd"))
    }, error = function(error_condition) {
      cat("Error in reading the dtm-model, file not found!")
      return()
    })
  useCoarseGrid <- FALSE
  # NORMALIZATION
  tryCatch(
    {
      cat("Normalizing height... ")
      sliVox_norm <- normalize_height(sliVox, dtm_z, na.rm = TRUE)
      cat("done!\n")
    }, error = function(error_condition) {
      useCoarseGrid <<- TRUE
      cat("-> problem with the normalisation, the min DTM model seems to be faulty...\n")
    })
  if(useCoarseGrid){
    # read in raster file 2, because first one produces some NAs at normalization
    cat("\"Roughly\" normalizing height - UNEXCEPTIONAL CASE!!! - (using the coarse 1 x 1 m dtm)...")
    tryCatch(
      {
        # reading raster file
        dtm_a <- raster(paste0(dirPath, groundPath, fileFinder, "_ground_rough.grd"))
      }, error = function(error_condition) {
        cat("Error in reading the rough dtm-model, file not found!")
        return()
      })

    tryCatch(
      {
        sliVox_norm <- normalize_height(sliVox, dtm_a, na.rm = TRUE) # need to save it in that intermediate object or it cannot unnormalize anymore
        cat("done!\n")
      }, error = function(error_condition) {
        cat("Error in creating the rough dtm model!")
        return()
      })
  }



  # WARNING: the variable sliVox always has real pointcloud and is original steep!
  # if you need the normalized point cloud, up until here it is sliVox_norm and just done five lines above.


  cat("Creating local output file slice_cluster_slope.laz containing", numClustCut, "clusters... ")
  # writeLAS(sliVox_norm, paste0(dbhPath, "slice_cluster.laz"))
  # global
  writeLAS(sliVox, paste0(dbhPath, "slice_cluster_slope.laz"))
  #local
  writeLAS(sliVox, paste0(dirPath, groundPath, fileFinder, "_clusterSlice_", bottomCut*100, "to", clipHeight*100, "_vox_slope.laz"))
  write.csv2(hilf, paste(dbhPath, "slice_cluster.csv", sep=""), row.names = F)
  sliVox <<- sliVox_norm
  v.env$sliVox <- sliVox_norm
  #points(tab.neu[tab.neu$cluster==0, ]$X, tab.neu[tab.neu$cluster==0, ]$Y, col=grey(0.8), pch=13, cex=0.3)
  #save.image(paste0(dbhPath, fileFinder, ".RData"))
  #write.csv(tab.neu, paste0(main_path, "Tabelle_nach_cluster_34.csv"))
  cat("done!\n")
  #}

  stop <- Sys.time()
  rm(sliVox_norm)
  gc()

  cat("Rough clustering is done.\n")
  print.difftime(round(stop - start, 1))

}







####################################################################################################
# 2 DIAMETER HUNT PER CLUSTER ######################################################################
####################################################################################################


diameterBeast_i <- function(clusterIndex, dbhPath, 
                            ipad = FALSE, fast = TRUE, allFiles = FALSE){
  
  # ERRORS INDUCED BY LEAVING EMPTY BRACKETS OVER WHOLE PARALLEL ROUTINE! ###
  suppressPackageStartupMessages(library("doParallel", character.only = T))
  suppressPackageStartupMessages(library("data.table", character.only = T))
  suppressPackageStartupMessages(library("ADPclust", character.only = T))
  suppressPackageStartupMessages(library("densityClust", character.only = T))
  #suppressPackageStartupMessages(library("dae", character.only = T))
  suppressPackageStartupMessages(library("plyr", character.only = T))
  suppressPackageStartupMessages(library("spatstat", character.only = T))
  suppressPackageStartupMessages(library("alphahull", character.only = T))
  suppressPackageStartupMessages(library("RANN", character.only = T))
  suppressPackageStartupMessages(library("flexclust", character.only = T))
  suppressPackageStartupMessages(library("sp", character.only = T))
  suppressPackageStartupMessages(library("matrixStats", character.only = T))
  #suppressPackageStartupMessages(library("Distance", character.only = T))
  suppressPackageStartupMessages(library("lmfor", character.only = T))
  suppressPackageStartupMessages(library("rgl", character.only = T))
  suppressPackageStartupMessages(library("conicfit", character.only = T))
  suppressPackageStartupMessages(library("MASS", character.only = T))
  suppressPackageStartupMessages(library("igraph", character.only = T))
  suppressPackageStartupMessages(library("geosphere", character.only = T))
  suppressPackageStartupMessages(library("pracma", character.only = T))
  suppressPackageStartupMessages(library("DescTools", character.only = T))
  suppressPackageStartupMessages(library("mgcv", character.only = T))
  suppressPackageStartupMessages(library("recexcavAAR", character.only = T))
  suppressPackageStartupMessages(library("raster", character.only = T))
  #suppressPackageStartupMessages(library("rlas", character.only = T))
  suppressPackageStartupMessages(library("lidR", character.only = T))
  suppressPackageStartupMessages(library("TreeLS", character.only = T))
  suppressPackageStartupMessages(library("dbscan", character.only = T))
  suppressPackageStartupMessages(library("rgl", character.only = T))
  suppressPackageStartupMessages(library("conicfit", character.only = T))
  #suppressPackageStartupMessages(library("VoxR", character.only = T))
  #suppressPackageStartupMessages(library("spatialEco", character.only = T))
  #suppressPackageStartupMessages(library("Rdistance", character.only = T))
  #suppressPackageStartupMessages(library("edci", character.only = T))
  #tryCatch(suppressPackageStartupMessages(library("edci")), 
  #         error = function(e){})
  
  #tab.neu <- tab.neu[tab.neu$cluster!=0, ] #alle weg die Noise sind
  u.grenzen.vec <- c(seq(1.0, 2.625, 0.125)) #Grenzen fuer BHD - Findung
  if(ipad){
    u.grenzen.vec <- c(seq(0.2, 1.825, 0.125)) #Grenzen fuer BHD - Findung #alt TLS/PLS: seq(1.0, 2.625, 0.125)
  }
  z.breite <- 0.15 #Breite der Schicht fuer circle/ell fit: 0.15
  #laenge.tab <- data.frame(id = Plot.ID.i, laenge = length(cluster.vec))
  
  path.output.cluster.end <- paste0(dbhPath, "fineCluster/")
  path.output.cluster.endgraph <- paste0(dbhPath, "graphSlice/")
  path.output.cluster.residuals <- paste0(dbhPath, "residuals/")
  
  if(allFiles){
    if(!dir.exists(path.output.cluster.endgraph)) dir.create(path.output.cluster.endgraph)
    if(!dir.exists(path.output.cluster.residuals)) dir.create(path.output.cluster.residuals)
  }
  
  angle_points <- function(xp, yp, xz, yz) {
    if(xp>xz&yp>yz|xp==xz&yp>yz|xp>xz&yp==yz){#1.Quadrant
      wink <- atan((xp-xz)/(yp-yz)) * 360/(2*pi)
    }
    if(xp>xz&yp<yz|xp==xz&yp<yz){#2.Quadrant
      wink <- atan((xp-xz)/(yp-yz)) * 360/(2*pi)+180
    }
    if(xp<xz&yp<yz){#3.Quadrant
      wink <- atan((xp-xz)/(yp-yz)) * 360/(2*pi)+180
    }
    if(xp<xz&yp>yz){#4.Quadrant
      wink <- atan((xp-xz)/(yp-yz)) * 360/(2*pi) + 360
    }
    if(xp<xz&yp==yz){#4.Quadrant - 270grad
      wink <- atan((xp-xz)/(yp-yz)) * 360/(2*pi) + 360
    }
    return(wink)
  }
  cat("~retr~")
  sliVox <- sliVox
  cat("~works~")
  cat(ifelse(exists("sliVox"), "existuje", "neex"))
  cat(ifelse(exists("sliVox"), 
             paste0("nPo", sliVox@header@PHB$`Number of point records`)))
  plot.clust2 <- filter_poi(sliVox, 'cluster' == clusterIndex)
  cat("~fí~")
  plot.clust2 <- data.frame("X" = plot.clust2$X, "Y" = plot.clust2$Y, "Z" = plot.clust2$Z, "cluster" = plot.clust2$cluster, "Intensity" = plot.clust2$Intensity)
  cat("~sec~")
  if(fast & nrow(plot.clust2) > 9000){
    # reduce stem to 9000 points to speed up diameter beast
    plot.clust2 <- plot.clust2[sample(nrow(plot.clust2), 9000, replace = F),]
  }
  
  cat("~tA~")
  #Herausgeben des richtigen Clusters in der jeweiligen Schicht
  #plot(plot.clust2)
  # plot3d(plot.clust2$X, plot.clust2$Y, plot.clust2$Z, col=plot.clust2$cluster, aspect=F)
  
  # hist(plot.clust2$Intensity)
  # mean(plot.clust2$Intensity)
  # median(plot.clust2$Intensity)
  #
  groesse <- (max(plot.clust2$X)-min(plot.clust2$X)) * (max(plot.clust2$Y)-min(plot.clust2$Y))
  groesse #m2
  
  cat("~tB~")
  if(groesse>=0.22){#0.22
    
    cat("~tC1~")
    hoehen <- seq(1.1, 2.5, 0.3)
    breite.hoehen <- 0.15
    test=4
    for(test in 1:length(hoehen)){
      plot.test <- plot.clust2[plot.clust2$Z>=hoehen[test]&plot.clust2$Z<=hoehen[test]+breite.hoehen, ]
      schoener.kreis <- data.frame("sK"=NA, "ksK"=NA)
      #plot(plot.test$X, plot.test$Y, asp=1)
      if(nrow(plot.test)!=0){
        par.ellipse <- EllipseDirectFit(cbind(plot.test$X, plot.test$Y))
        geom.ellipse <- as.vector(AtoG(par.ellipse)$ParG)
        
        cat("~tCloop~")
        n.vertices <- 400
        
        ellipse.vert <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3], b=geom.ellipse[4],
                                         angle=180/pi*geom.ellipse[5], steps=n.vertices)
        
        q <- geom.ellipse[4]/geom.ellipse[3]
        if(q>1){q <- geom.ellipse[3]/geom.ellipse[4]}
        x.ellipse <- geom.ellipse[1]
        y.ellipse <- geom.ellipse[2]
        dbh.ellipse <- round(sqrt(((2*geom.ellipse[4])^2+(2*geom.ellipse[3])^2)/2) * 100, 2)
        #points(ellipse.vert, type="l", lwd=3, lty=2, col=4)
        
        # #Puffer fuer Ellipse
        ellipse.vert.gross <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3]*1.1, b=geom.ellipse[4]*1.15,
                                               angle=180/pi*geom.ellipse[5], steps=n.vertices)
        
        ellipse.vert.klein <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3]*0.85, b=geom.ellipse[4]*0.85,
                                               angle=180/pi*geom.ellipse[5], steps=n.vertices)
        
        klein_ell <- point.in.polygon(plot.test$X, plot.test$Y, ellipse.vert.klein[, 1], ellipse.vert.klein[, 2]) #wieviele Punkte sind in der kleinen Ell
        gross_ell <- point.in.polygon(plot.test$X, plot.test$Y, ellipse.vert.gross[, 1], ellipse.vert.gross[, 2]) #wieviele in der grossen
        
        #points(ellipse.vert.klein, type="l", lwd=3, lty=2, col=2)
        #points(ellipse.vert.gross, type="l", lwd=3, lty=2, col=3)
        anteil_punkte_ell <- sum(klein_ell)/sum(gross_ell)
        anteil_punkte_segment <- (sum(gross_ell)-sum(klein_ell))/nrow(plot.test)
        
        if(anteil_punkte_ell<=0.45&anteil_punkte_segment>0.42&q>0.75){
          schoener.kreis$sK <- 1
          schoener.kreis$ksK <- 0
          schoener.kreis$x <- x.ellipse
          schoener.kreis$y <- y.ellipse
          schoener.kreis$BHD <- dbh.ellipse
          schoener.kreis$q <- q
        }else{
          schoener.kreis$sK <- 0
          schoener.kreis$ksK <- 1
          schoener.kreis$x <- x.ellipse
          schoener.kreis$y <- y.ellipse
          schoener.kreis$BHD <- dbh.ellipse
          schoener.kreis$q <- q
        }
      }else{
        schoener.kreis$sK <- 0
        schoener.kreis$ksK <- 1
        schoener.kreis$x <- 1
        schoener.kreis$y <- 1
        schoener.kreis$BHD <- 1
        schoener.kreis$q <- 0.2
      }
      
      
      if(test==1){
        schoener.kreis.out <- schoener.kreis
      }else{
        schoener.kreis.out <- rbind(schoener.kreis.out, schoener.kreis)
      }
      
    }
    
    cat("~NEVER~")
    schoener.kreis.out
    mean(schoener.kreis.out$q)
    sd(schoener.kreis.out$BHD)
    if(sum(schoener.kreis.out$sK)<3&sd(schoener.kreis.out$BHD)>1.5&mean(schoener.kreis.out$q)<0.9|groesse>1.1|groesse>0.4&sum(schoener.kreis.out$sK)<3){
      plot.clust2 <- filter_poi(sliVox, cluster == clusterIndex) #Herausgeben des richtigen Clusters in der jeweiligen Schicht
      plot.clust2 <- data.frame("X" = plot.clust2$X, "Y" = plot.clust2$Y, "Z" = plot.clust2$Z, "cluster" = plot.clust2$cluster, "Intensity" = plot.clust2$Intensity)
      
      zalt <- c(1, 3)
      zneu <- c(1.8, 2.2)
      model <- lm(zneu~zalt)
      plot.clust2$Zalt<- plot.clust2$Z
      plot.clust2$Z <- predict(model, newdata=data.frame("zalt"=plot.clust2$Zalt))
      #2 in 314 braechte 22; 8_1 auch Problem
      res <- optics(plot.clust2[, 1:3], eps = 0.025,  minPts = 20) #eps = 0.025,  minPts = 21
      res
      #Cluster herausgeben
      res <- extractDBSCAN(res, eps_cl = .02) #0.02
      res
      plot.clust2$cluster <- res$cluster
      plot.clust2 <- plot.clust2[plot.clust2$cluster!=0, ]
      #plot3d(plot.clust2$X, plot.clust2$Y, plot.clust2$Z, col=plot.clust2$cluster, aspect = F)
      plot.clust2$Z <- plot.clust2$Zalt
      
      hilf.tab <- as.data.frame(table(res$cluster))
      colnames(hilf.tab) <- c("cluster", "anz")
      hilf.tab <- hilf.tab[hilf.tab$cluster!=0, ]
      hilf.tab <- hilf.tab[order(-hilf.tab$anz), ]
      
      if(nrow(hilf.tab)>=10){
        zahl_cluster <- hilf.tab[1:10, ]$cluster
        plot.clust2 <- plot.clust2[plot.clust2$cluster%in%zahl_cluster, ]
      }else{
        zahl_cluster <- hilf.tab[1:nrow(hilf.tab), ]$cluster
        plot.clust2 <- plot.clust2[plot.clust2$cluster%in%zahl_cluster, ]
      }
      clu2 <- 1
      for(clu2 in 1:length(unique(plot.clust2$cluster))){ #bei den 3 mit den meisten Punkten schauen, ob die Baeume sind mit vertikaler ausdehnung
        clu.i <- plot.clust2[plot.clust2$cluster==unique(plot.clust2$cluster)[clu2], ]
        ausdehnung_z <- max(clu.i$Z)-min(clu.i$Z)
        if(ausdehnung_z>1.2){
          ausdehnung_out <- clu.i
        }else{
          ausdehnung_out <- data.frame()
        }
        if(clu2==1){
          plot.clu <- ausdehnung_out
        }else{
          plot.clu <- rbind(plot.clu, ausdehnung_out)
        }
      }
      
      plot.clust2 <- plot.clu
      #table(plot.clust2$cluster)
      
      #plot3d(plot.clust2$X, plot.clust2$Y, plot.clust2$Z, col=plot.clust2$cluster, aspect=F)
      
      if(nrow(plot.clust2)<500&groesse<0.6){
        plot.clust2 <- filter_poi(sliVox, cluster == clusterIndex) #Herausgeben des richtigen Clusters in der jeweiligen Schicht
        plot.clust2 <- data.frame("X" = plot.clust2$X, "Y" = plot.clust2$Y, "Z" = plot.clust2$Z, "cluster" = plot.clust2$cluster, "Intensity" = plot.clust2$Intensity)
      }
    }
    
  }
  else{
    zalt <- c(1, 3)
    zneu <- c(1.8, 2.2)
    model <- lm(zneu~zalt)
    plot.clust2$Zalt<- plot.clust2$Z
    plot.clust2$Z <- predict(model, newdata=data.frame("zalt"=plot.clust2$Zalt))
    
    res <- optics(plot.clust2[, 1:3], eps = 0.025,  minPts = 18) #eps = 0.025,  minPts = 18
    res
    #Cluster herausgeben
    res <- extractDBSCAN(res, eps_cl = .023) #0.02
    res
    plot.clust2$cluster <- res$cluster
    plot.clust2 <- plot.clust2[plot.clust2$cluster!=0, ]
    #plot3d(plot.clust2$X, plot.clust2$Y, plot.clust2$Z, col=plot.clust2$cluster, aspect = F)
    plot.clust2$Z <- plot.clust2$Zalt
    
    hilf.tab <- as.data.frame(table(res$cluster))
    colnames(hilf.tab) <- c("cluster", "anz")
    hilf.tab <- hilf.tab[hilf.tab$cluster!=0, ]
    hilf.tab <- hilf.tab[order(-hilf.tab$anz), ]
    
    if(nrow(hilf.tab)>=10){
      zahl_cluster <- hilf.tab[1:10, ]$cluster
      plot.clust2 <- plot.clust2[plot.clust2$cluster%in%zahl_cluster, ]
    }else{
      zahl_cluster <- hilf.tab[1:nrow(hilf.tab), ]$cluster
      plot.clust2 <- plot.clust2[plot.clust2$cluster%in%zahl_cluster, ]
    }
    clu2 <- 1
    for(clu2 in 1:length(unique(plot.clust2$cluster))){ #bei den 3 mit den meisten Punkten schauen, ob die Baeume sind mit vertikaler ausdehnung
      clu.i <- plot.clust2[plot.clust2$cluster==unique(plot.clust2$cluster)[clu2], ]
      ausdehnung_z <- max(clu.i$Z)-min(clu.i$Z)
      if(ausdehnung_z>1.2){
        ausdehnung_out <- clu.i
      }else{
        ausdehnung_out <- data.frame()
      }
      if(clu2==1){
        plot.clu <- ausdehnung_out
      }else{
        plot.clu <- rbind(plot.clu, ausdehnung_out)
      }
    }
    
    plot.clust2 <- plot.clu
    #table(plot.clust2$cluster)
    
    #plot3d(plot.clust2$X, plot.clust2$Y, plot.clust2$Z, col=plot.clust2$cluster, aspect=F)
    
    if(nrow(plot.clust2)<500){
      plot.clust2 <- filter_poi(sliVox, cluster == clusterIndex) #Herausgeben des richtigen Clusters in der jeweiligen Schicht
      plot.clust2 <- data.frame("X" = plot.clust2$X, "Y" = plot.clust2$Y, "Z" = plot.clust2$Z, "cluster" = plot.clust2$cluster, "Intensity" = plot.clust2$Intensity)
    }
  }
  
  cluster2 <- unique(plot.clust2$cluster)
  #plot3d(plot.clust2$X, plot.clust2$Y, plot.clust2$Z, col=plot.clust2$cluster, aspect=F)
  
  clust2=1
  j=1
  
  if(nrow(plot.clust2)==0){
    return()
  }
  
  for(clust2 in 1:length(cluster2)){
    for(j in 1:length(u.grenzen.vec)){
      plot.i <- plot.clust2[plot.clust2$cluster==cluster2[clust2]&plot.clust2$Z>=u.grenzen.vec[j]&plot.clust2$Z<=u.grenzen.vec[j]+z.breite, ] #Herausgeben des richtigen Clusters in der jeweiligen Schicht
      if(nrow(plot.i)>5){
        par.ellipse <- EllipseDirectFit(cbind(plot.i$X, plot.i$Y))
        geom.ellipse <- as.vector(AtoG(par.ellipse)$ParG)
        q <- geom.ellipse[4]/geom.ellipse[3]
        if(q>1){q <- geom.ellipse[3]/geom.ellipse[4]}
      }else{
        q <- NA
      }
      
      abflachung <- data.frame("cluster"=as.numeric(cluster2[clust2]), "q"=q)
      if(j==1){
        abflachung1 <- abflachung
      }else{
        abflachung1 <- rbind(abflachung1, abflachung)
      }
    }
    med <- median(abflachung1$q, na.rm = TRUE)
    if(!is.na(med)&med<0.65&nrow(plot.clust2)>1000){
      j=2
      for(j in 1:length(u.grenzen.vec)){
        plot.i <- plot.clust2[plot.clust2$cluster==cluster2[clust2]&plot.clust2$Z>=u.grenzen.vec[j]&plot.clust2$Z<=u.grenzen.vec[j]+z.breite, ] #Herausgeben des richtigen Clusters in der jeweiligen Schicht
        #plot(plot.i$X, plot.i$Y, asp=1)
        if(!nrow(plot.i)==0){
          res <- optics(plot.i[, 1:2], eps = 0.03,  minPts = 30) #eps = 0.03,  minPts = 30
          res
          
          #Cluster herausgeben
          res <- extractDBSCAN(res, eps_cl = .022) #0.022
          res
          
          plot.i$cluster2 <- res$cluster
          hilf.tab <- as.data.frame(table(res$cluster))
          colnames(hilf.tab) <- c("cluster", "anz")
          hilf.tab <- hilf.tab[hilf.tab$cluster!=0, ]
          
          ###
          # new 500 says CG 2022-08-22
          #zahl_cluster <- hilf.tab[which(hilf.tab$anz>=100), ]$cluster
          zahl_cluster <- hilf.tab[which(hilf.tab$anz>=500), ]$cluster
          plot.i <- plot.i[plot.i$cluster2%in%zahl_cluster, ]
          #plot(plot.i$X, plot.i$Y, asp=1, col=plot.i$cluster2)
          
          if(length(zahl_cluster)>=2){
            zwei <- 1
          }else{
            zwei <- 0
          }
          
        }else{
          zwei <- 0
          plot.i$cluster2 <- NA
        }
        
        plot.zwei <- plot.i
        plot.zwei$ho <- u.grenzen.vec[j]
        plot.zwei$zwei <- zwei
        
        if(j==1){
          plot.zwei.out <- plot.zwei
          zwei.out <- zwei
        }else{
          plot.zwei.out <- rbind(plot.zwei.out, plot.zwei)
          zwei.out <- c(zwei.out, zwei)
        }
      }
      if(sum(zwei.out)>=6){
        zwei_baum <- plot.zwei.out[plot.zwei.out$zwei==1, ]
        
        #plot3d(zwei_baum$X, zwei_baum$Y, zwei_baum$Z, asp=F)
        
        zalt <- c(1, 3)
        zneu <- c(1.8, 2.2)
        model <- lm(zneu~zalt)
        zwei_baum$Zalt<- zwei_baum$Z
        zwei_baum$Z <- predict(model, newdata=data.frame("zalt"=zwei_baum$Zalt))
        
        res <- optics(zwei_baum[, 1:3], eps = 0.025,  minPts = 22) #eps = 0.025,  minPts = 18
        res
        #Cluster herausgeben
        res <- extractDBSCAN(res, eps_cl = .020) #0.02
        res
        zwei_baum$cluster2 <- res$cluster
        zwei_baum <- zwei_baum[zwei_baum$cluster2!=0, ]
        #plot3d(zwei_baum$X, zwei_baum$Y, zwei_baum$Z, col=zwei_baum$cluster, aspect = F)
        
        zwei_baum$Z <- zwei_baum$Zalt
        
        hilf.tab <- as.data.frame(table(res$cluster))
        colnames(hilf.tab) <- c("cluster", "anz")
        hilf.tab <- hilf.tab[hilf.tab$cluster!=0, ]
        hilf.tab <- hilf.tab[order(-hilf.tab$anz), ]
        
        if(nrow(hilf.tab)>=10){
          zahl_cluster <- hilf.tab[1:10, ]$cluster
          zwei_baum <- zwei_baum[zwei_baum$cluster2%in%zahl_cluster, ]
        }else{
          zahl_cluster <- hilf.tab[1:nrow(hilf.tab), ]$cluster
          zwei_baum <- zwei_baum[zwei_baum$cluster2%in%zahl_cluster, ]
        }
        clu2 <- 2
        for(clu2 in 1:length(unique(zwei_baum$cluster2))){ #bei den 3 mit den meisten Punkten schauen, ob die Baeume sind mit vertikaler ausdehnung
          clu.i <- zwei_baum[zwei_baum$cluster2==unique(zwei_baum$cluster2)[clu2], ]
          ausdehnung_z <- max(clu.i$Z)-min(clu.i$Z)
          if(ausdehnung_z>0.7){
            ausdehnung_out <- clu.i
          }else{
            ausdehnung_out <- data.frame()
          }
          if(clu2==1){
            plot.clu <- ausdehnung_out
          }else{
            plot.clu <- rbind(plot.clu, ausdehnung_out)
          }
        }
        
        plot.clust_end <- plot.clu
        plot.clust_end <- plot.clust_end[, c("X", "Y", "Z", "Intensity", "cluster", "cluster2")]
        if(nrow(plot.clust_end)<300){
          plot.clust_end <- plot.clust2[plot.clust2$cluster==cluster2[clust2], ]
          plot.clust_end$cluster2 <- 333 #0
          plot.clust_end <- plot.clust_end[, c("X", "Y", "Z", "Intensity", "cluster", "cluster2")]
        }
      }else{
        plot.clust_end <- plot.clust2[plot.clust2$cluster==cluster2[clust2], ]
        plot.clust_end$cluster2 <- 333 #0
        plot.clust_end <- plot.clust_end[, c("X", "Y", "Z", "Intensity", "cluster", "cluster2")]
      }
    }else{
      plot.clust_end <- plot.clust2[plot.clust2$cluster==cluster2[clust2], ]
      plot.clust_end$cluster2 <- 333 #0
      plot.clust_end <- plot.clust_end[, c("X", "Y", "Z", "Intensity", "cluster", "cluster2")]
    }
    
    
    if(clust2==1){
      plot.clust_final <- plot.clust_end
    }else{
      plot.clust_final <- rbind(plot.clust_final, plot.clust_end)
    }
  }
  
  
  plot.clust2 <- plot.clust_final
  (cluster2 <- as.numeric(unique(plot.clust2$cluster)))
  #(cluster3 <- unique(plot.clust2$cluster2))
  # plot3d(plot.clust2$X, plot.clust2$Y, plot.clust2$Z, col=plot.clust2$cluster2, aspect=F)
  # hist(plot.clust2$Intensity)
  # mean(plot.clust2$Intensity)
  # median(plot.clust2$Intensity)
  
  (quantil_inten <- quantile(plot.clust2$Intensity, prob=0.8)) #muss ueber 10000 sein
  
  if(ipad){
    quantil_inten <- 8000
  }
  
  plot.clust2$cluster <- as.numeric(plot.clust2$cluster)
  plot.clust2$cluster2 <- as.numeric(plot.clust2$cluster2)
  
  class(plot.clust2$cluster2)
  class(cluster2[clust2])
  clust2=1
  clust3=1
  j=1
  # says CG 2022-08-22 new but not accepted
  # if(nrow(plot.clust2)!=0&quantil_inten>=7900){ #mind 7900 muss die intensitaet quantile sein
  if(nrow(plot.clust2)!=0&quantil_inten>=6000){ #vormals 7900 #jetzt 6000
    
    for(clust2 in 1:length(cluster2)){
      nr <- cluster2[clust2]
      plot.clust3 <- plot.clust2[plot.clust2$cluster==nr, ] #Herausgeben des richtigen Clusters in der jeweiligen Schicht
      (cluster3 <- unique(plot.clust3$cluster2))
      for(clust3 in 1:length(cluster3)){
        for(j in 1:length(u.grenzen.vec)){
          plot.i <- plot.clust3[plot.clust3$cluster2==cluster3[clust3]&plot.clust3$Z>=u.grenzen.vec[j]&plot.clust3$Z<=u.grenzen.vec[j]+z.breite, ] #Herausgeben des richtigen Clusters in der jeweiligen Schicht
          if(nrow(plot.i)>4){
            #plot(plot.i$X, plot.i$Y, asp=1)
            nrow(plot.i)
            #plot(plot.i$X, plot.i$Y, asp=1) #Zeichnen der Punkte
            
            par.ellipse <- EllipseDirectFit(cbind(plot.i$X, plot.i$Y))
            geom.ellipse <- as.vector(AtoG(par.ellipse)$ParG)
            q <- geom.ellipse[4]/geom.ellipse[3]
            if(q>1){q <- geom.ellipse[3]/geom.ellipse[4]}
            n.vertices <- 400
            ellipse.gross <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3]*1.1, b=geom.ellipse[4]*1.1,
                                              angle=180/pi*geom.ellipse[5], steps=n.vertices)
            ellipse.klein <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3]*0.9, b=geom.ellipse[4]*0.9,
                                              angle=180/pi*geom.ellipse[5], steps=n.vertices)
            
            p_klein_ell <- point.in.polygon(plot.i$X, plot.i$Y, ellipse.klein[, 1], ellipse.klein[, 2])
            p_gross_ell <- point.in.polygon(plot.i$X, plot.i$Y, ellipse.gross[, 1], ellipse.gross[, 2])
            
            anteil_punk_ell <- sum(p_klein_ell)/sum(p_gross_ell)
            #points(ellipse.gross, type="l", lwd=3, lty=2, col=2)
            #points(ellipse.klein, type="l", lwd=3, lty=2, col=3)
            
            ell.segment <- plot.i[p_gross_ell==1&p_klein_ell!=1, ]
            if(nrow(ell.segment)!=0){
              ell.segment$winkel <- 1
              neuew=1
              for(neuew in 1:nrow(ell.segment)){
                ell.segment[neuew, "winkel"] <- angle_points(ell.segment[neuew, "X"], ell.segment[neuew, "Y"],
                                                             geom.ellipse[1], geom.ellipse[2])
              }
              winkel <- seq(0, 350, 10)
              breite <- 10
              wink=1
              for(wink in 1:length(winkel)){
                untere_grenz <- winkel[wink]
                obere_grenz <- winkel[wink]+breite
                sub.i <- ell.segment[ell.segment$winkel>=untere_grenz&ell.segment$winkel<obere_grenz, ]
                if(!nrow(sub.i)==0){
                  gefuellt <- 1
                }else{
                  gefuellt <- 0
                }
                if(wink==1){
                  gefuellt.out <- gefuellt
                }else{
                  gefuellt.out <- c(gefuellt.out, gefuellt)}
              }
              prozent_fuellung <- sum(gefuellt.out)/length(gefuellt.out)*100
            }else{
              prozent_fuellung <- 0
            }
            
            
            if(nrow(plot.i)>50){ #50
              if(nrow(plot.i)>1300&!q>0.88&!prozent_fuellung>85|nrow(plot.i)>1300&anteil_punk_ell>0.55){
                #Ordering points to identify the clustering structure
                res <- optics(plot.i[, 1:2], eps = 0.03,  minPts = 30) #eps = 0.03, minPts=40
                res
                
                #Cluster herausgeben
                res <- extractDBSCAN(res, eps_cl = .025) #0.2 #0.1 auch gut/besser 0.04-0.08 #0.025
                res
                
                plot.i$cluster2 <- res$cluster
                #plot(plot.i$X, plot.i$Y, asp=1, col=plot.i$cluster2)
                
                hilf.tab <- as.data.frame(table(res$cluster))
                colnames(hilf.tab) <- c("cluster", "anz")
                hilf.tab <- hilf.tab[hilf.tab$cluster!=0, ]
                
                ###
                zahl_cluster <- hilf.tab[which(hilf.tab$anz==max(hilf.tab$anz)), ]$cluster[1]
                
                #plot(plot.i$X, plot.i$Y, col=plot.i$cluster2, asp=1, cex=0.5) #Zeichnen der geclusterten Punkte
                
                if(nrow(hilf.tab)!=0){
                  
                  #Ellipsen Fit
                  plot.ell <- plot.i[plot.i$cluster2!=0, ]
                  #plot.ell <- plot.i
                  par.ellipse <- EllipseDirectFit(cbind(plot.ell$X, plot.ell$Y))
                  geom.ellipse <- as.vector(AtoG(par.ellipse)$ParG)
                  
                  n.vertices <- 400
                  
                  ellipse.vert <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3], b=geom.ellipse[4],
                                                   angle=180/pi*geom.ellipse[5], steps=n.vertices)
                  
                  q <- geom.ellipse[4]/geom.ellipse[3]
                  if(q>1){q <- geom.ellipse[3]/geom.ellipse[4]}
                  x.ellipse <- geom.ellipse[1]
                  y.ellipse <- geom.ellipse[2]
                  dbh.ellipse <- round(sqrt(((2*geom.ellipse[4])^2+(2*geom.ellipse[3])^2)/2) * 100, 2)
                  #points(ellipse.vert, type="l", lwd=3, lty=2, col=4)
                  
                  # #Puffer fuer Ellipse
                  ellipse.vert.gross <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3]*1.1, b=geom.ellipse[4]*1.1,
                                                         angle=180/pi*geom.ellipse[5], steps=n.vertices)
                  
                  ellipse.vert.klein <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3]*0.85, b=geom.ellipse[4]*0.85,
                                                         angle=180/pi*geom.ellipse[5], steps=n.vertices)
                  
                  klein_ell <- point.in.polygon(plot.ell$X, plot.ell$Y, ellipse.vert.klein[, 1], ellipse.vert.klein[, 2]) #wieviele Punkte sind in der kleinen Ell
                  gross_ell <- point.in.polygon(plot.ell$X, plot.ell$Y, ellipse.vert.gross[, 1], ellipse.vert.gross[, 2]) #wieviele in der grossen
                  
                  #points(ellipse.vert.klein, type="l", lwd=3, lty=2, col=2)
                  #points(ellipse.vert.gross, type="l", lwd=3, lty=2, col=3)
                  anteil_punkte_ell <- sum(klein_ell)/sum(gross_ell)
                  
                  if(anteil_punkte_ell<=0.45&q>0.55){
                    plot.i.ell.out <- plot.i
                    plot.i.ell.out <- plot.i.ell.out[plot.i.ell.out$cluster2!=zahl_cluster, ]
                    
                    set.seed(15)
                    if( (class(try(
                      find.clust.circ <- circMclust( datax=plot.i$X,
                                                     datay=plot.i$Y,
                                                     method="const",
                                                     # nx=1, ny=1, nr=1, # nx=10, ny=10, nr=25,
                                                     nx=8, ny=8, nr=4, # nx=10, ny=10, nr=5,
                                                     minsr=0.01, maxsr=0.5,
                                                     nc=1, # nc=1,
                                                     # minsd=0.1, maxsd=0.5,
                                                     bw=0.01 # bw=0.05
                      )
                      
                    ))=="try-error")==FALSE ){
                      try.error <- FALSE
                    }else{
                      try.error <- TRUE
                    }
                    
                    
                    if(is.null(find.clust.circ)&try.error==FALSE){ #Problem beim Circle Fit abfangen
                      dbh.circle <- NA
                      x.circle <- NA
                      y.circle <- NA
                      anteil_punkte_circ <- NA
                    }
                    
                    plot.i <- plot.i[plot.i$cluster2==zahl_cluster, ]
                    
                    if(!is.null(find.clust.circ)&try.error==F){ #Problem beim Circle Fit abfangen
                      if(colnames(find.clust.circ)[4]=="value"){ #Problem beim Circle Fit abfangen
                        if(is.na(as.numeric(bestMclust(find.clust.circ, 1))[1])==FALSE){ #Problem beim Circle Fit abfangen
                          par.circle.1 <- as.numeric(bestMclust(find.clust.circ, 1))
                          dbh.circle <- round(2 * par.circle.1[3] * 100, 2)
                          x.circle <- par.circle.1[1]
                          y.circle <- par.circle.1[2]
                          
                          # points(par.circle.1[1], par.circle.1[2], cex=1.5, col=6, pch=13)
                          # circle(par.circle.1[1], par.circle.1[2], par.circle.1[3], border=6, lwd=2)
                          
                          circle.vert.gross <- calculateCircle(par.circle.1[1], par.circle.1[2], par.circle.1[3]*1.2, steps=n.vertices)
                          circle.vert.klein <- calculateCircle(par.circle.1[1], par.circle.1[2], par.circle.1[3]*0.8, steps=n.vertices)
                          
                          kleine <- point.in.polygon(plot.i.ell.out$X, plot.i.ell.out$Y, circle.vert.klein[, 1], circle.vert.klein[, 2])
                          grosse <- point.in.polygon(plot.i.ell.out$X, plot.i.ell.out$Y, circle.vert.gross[, 1], circle.vert.gross[, 2])
                          
                          plot.i.ell.out <- plot.i.ell.out[grosse==1&kleine!=1, ]
                          #points(plot.i.ell.out$X, plot.i.ell.out$Y, col=2)
                          
                          plot.i <- rbind(plot.i, plot.i.ell.out)
                          schoener.kreis <- "ja"
                        }else{
                          plot.i <- plot.i[plot.i$cluster2==zahl_cluster, ]
                          schoener.kreis <- "ja"}
                      }else{
                        plot.i <- plot.i[plot.i$cluster2==zahl_cluster, ]
                        schoener.kreis <- "ja"}
                    }else{
                      plot.i <- plot.i[plot.i$cluster2==zahl_cluster, ]
                      schoener.kreis <- "ja"}
                    
                    
                  }else{
                    if(nrow(plot.i[plot.i$cluster2==zahl_cluster, ])>1400){ #1200
                      plot.i.new.cluster <- plot.i[plot.i$cluster2==zahl_cluster, ]
                      
                      #Ordering points to identify the clustering structure
                      res <- optics(plot.i.new.cluster[, 1:2], eps = 0.03,  minPts = 100) #eps = 0.2, minPts=50
                      res
                      
                      #Cluster herausgeben
                      res <- extractDBSCAN(res, eps_cl = .02) #0.2 #0.1 auch gut/besser 0.04-0.08 #0.025
                      res
                      
                      plot.i.new.cluster$cluster2 <- res$cluster
                      #plot(plot.i.new.cluster$X, plot.i.new.cluster$Y, asp=1, col=plot.i.new.cluster$cluster2)
                      
                      hilf.tab <- as.data.frame(table(res$cluster))
                      colnames(hilf.tab) <- c("cluster", "anz")
                      hilf.tab <- hilf.tab[hilf.tab$cluster!=0, ]
                      
                      ###
                      zahl_cluster <- hilf.tab[which(hilf.tab$anz==max(hilf.tab$anz)), ]$cluster[1]
                      
                      plot.i <- plot.i.new.cluster[plot.i.new.cluster$cluster2==zahl_cluster, ]
                      schoener.kreis <- "nein"
                    }else{
                      plot.i.new.cluster <- plot.i[plot.i$cluster2==zahl_cluster, ]
                      
                      #Ordering points to identify the clustering structure
                      res <- optics(plot.i.new.cluster[, 1:2], eps = 0.03,  minPts = 40) #eps = 0.2, minPts=50
                      res
                      
                      #Cluster herausgeben
                      res <- extractDBSCAN(res, eps_cl = .025) #0.2 #0.1 auch gut/besser 0.04-0.08 #0.025
                      res
                      
                      plot.i.new.cluster$cluster2 <- res$cluster
                      #plot(plot.i.new.cluster$X, plot.i.new.cluster$Y, asp=1, col=plot.i.new.cluster$cluster2)
                      
                      hilf.tab <- as.data.frame(table(res$cluster))
                      colnames(hilf.tab) <- c("cluster", "anz")
                      hilf.tab <- hilf.tab[hilf.tab$cluster!=0, ]
                      
                      ###
                      zahl_cluster <- hilf.tab[which(hilf.tab$anz==max(hilf.tab$anz)), ]$cluster[1]
                      
                      plot.i <- plot.i.new.cluster[plot.i.new.cluster$cluster2==zahl_cluster, ]
                      schoener.kreis <- "nein"
                    }
                    
                    
                    
                  }
                  
                  
                  if(nrow(plot.i)>120){ #nur fitten wenn mehr als 5 Punkte da sind
                    
                    anzahl.punkte <- nrow(plot.i)
                    
                    plot.j <- plot.i
                    # ?circMclust
                    if(allFiles){
                      png(paste(path.output.cluster.endgraph, "cluster_", clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], "_", "hoehe_", u.grenzen.vec[j] , ".png", sep=""),
                          type = "cairo", width = 1000, height = 1000)
                      
                      par(mfrow=c(1, 1))
                    }
                    
                    #plot(plot.j$X, plot.j$Y, asp=1) #Zeichnen der Punkte
                    
                    #Endgueltige Ellipse
                    par.ellipse <- EllipseDirectFit(cbind(plot.j$X, plot.j$Y))
                    geom.ellipse <- as.vector(AtoG(par.ellipse)$ParG)
                    n.vertices <- 400
                    ellipse.vert <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3], b=geom.ellipse[4],
                                                     angle=180/pi*geom.ellipse[5], steps=n.vertices)
                    q <- geom.ellipse[4]/geom.ellipse[3]
                    if(q>1){q <- geom.ellipse[3]/geom.ellipse[4]}
                    x.ellipse <- geom.ellipse[1]
                    y.ellipse <- geom.ellipse[2]
                    dbh.ellipse <- sqrt(((2*geom.ellipse[4])^2+(2*geom.ellipse[3])^2)/2) * 100
                    
                    
                    #Endgueltiger Kreis
                    #Circle Fit
                    set.seed(15)
                    if( (class(try(
                      if(fast){
                        find.clust.circ <- circMclust( datax=plot.j$X,
                                                       datay=plot.j$Y,
                                                       method="const",
                                                       # nx=1, ny=1, nr=1, # nx=10, ny=10, nr=25,
                                                       nx=8, ny=8, nr=4, #nx=25, ny=25, nr=5,
                                                       minsr=0.01, maxsr=0.5,
                                                       nc=1, # nc=1,
                                                       # minsd=0.1, maxsd=0.5,
                                                       bw=0.01 # bw=0.05
                                                       
                        )
                      } else {
                        # extended nx, ny and nr, takes way longer! would need more extensive testing, but before other efficiency measures
                        find.clust.circ <- circMclust( datax=plot.j$X,
                                                       datay=plot.j$Y,
                                                       method="const",
                                                       # nx=1, ny=1, nr=1, # nx=10, ny=10, nr=25,
                                                       nx=10, ny=10, nr=5, # nx=25, ny=25, nr=5, MAX DETAIL
                                                       minsr=0.01, maxsr=0.5,
                                                       nc=1, # nc=1,
                                                       # minsd=0.1, maxsd=0.5,
                                                       bw=0.01 # bw=0.05
                        )
                      }
                    ))=="try-error")==FALSE ){
                      try.error <- FALSE
                    }else{
                      try.error <- TRUE
                    }
                    
                    
                    if(is.null(find.clust.circ)&try.error==FALSE){ #Problem beim Circle Fit abfangen
                      dbh.circle <- NA
                      x.circle <- NA
                      y.circle <- NA
                      anteil_punkte_circ <- NA
                    }
                    
                    if(!is.null(find.clust.circ)&try.error==F){ #Problem beim Circle Fit abfangen
                      if(colnames(find.clust.circ)[4]=="value"){ #Problem beim Circle Fit abfangen
                        if(is.na(as.numeric(bestMclust(find.clust.circ, 1))[1])==FALSE){ #Problem beim Circle Fit abfangen
                          par.circle.1 <- as.numeric(bestMclust(find.clust.circ, 1))
                          dbh.circle <- round(2 * par.circle.1[3] * 100, 6)
                          x.circle <- par.circle.1[1]
                          y.circle <- par.circle.1[2]
                          
                          #points(par.circle.1[1], par.circle.1[2], cex=1.5, col=6, pch=13)
                          #circle(par.circle.1[1], par.circle.1[2], par.circle.1[3], border=6, lwd=2)
                          
                          #Zuschlag fuer Ausschneiden
                          if(dbh.circle<30){zuschlag=0.06}else{zuschlag=0.09}
                          
                          circle.vert.gross <- calculateCircle(par.circle.1[1], par.circle.1[2], par.circle.1[3]+zuschlag, steps=n.vertices)
                          circle.vert.klein <- calculateCircle(par.circle.1[1], par.circle.1[2], par.circle.1[3]*0.8, steps=n.vertices)
                          
                          klein_circ <- point.in.polygon(plot.j$X, plot.j$Y, circle.vert.klein[, 1], circle.vert.klein[, 2])
                          gross_circ <- point.in.polygon(plot.j$X, plot.j$Y, circle.vert.gross[, 1], circle.vert.gross[, 2])
                          
                          # points(circle.vert.klein, type="l", lwd=3, lty=2, col=2)
                          # points(circle.vert.gross, type="l", lwd=3, lty=2, col=3)
                          
                          anteil_punkte_circ <- sum(klein_circ)/sum(gross_circ)
                          
                          kreis.segment <- plot.j[gross_circ==1&klein_circ!=1, ]
                          #points(kreis.segment$X, kreis.segment$Y, asp=1, col=2)
                          kreis.segment$winkel <- 1
                          neuew=1
                          for(neuew in 1:nrow(kreis.segment)){
                            kreis.segment[neuew, "winkel"] <- angle_points(kreis.segment[neuew, "X"], kreis.segment[neuew, "Y"],
                                                                           x.circle, y.circle)
                          }
                          
                          #hist(kreis.segment$winkel, breaks=seq(0, 350, 10))
                          winkel <- seq(0, 350, 10)
                          breite <- 10
                          wink=1
                          
                          for(wink in 1:length(winkel)){
                            untere_grenz <- winkel[wink]
                            obere_grenz <- winkel[wink]+breite
                            sub.i <- kreis.segment[kreis.segment$winkel>=untere_grenz&kreis.segment$winkel<obere_grenz, ]
                            if(!nrow(sub.i)==0){
                              gefuellt <- 1
                            }else{
                              gefuellt <- 0
                            }
                            if(wink==1){
                              gefuellt.out <- gefuellt
                            }else{
                              gefuellt.out <- c(gefuellt.out, gefuellt)}
                          }
                          
                          prozent_fuellung <- sum(gefuellt.out)/length(gefuellt.out)*100
                          
                          
                        }else{
                          dbh.circle <- NA
                          x.circle <- NA
                          y.circle <- NA
                          anteil_punkte_circ <- NA
                          prozent_fuellung <- NA
                        }
                        
                      }else{
                        dbh.circle <- NA
                        x.circle <- NA
                        y.circle <- NA
                        anteil_punkte_circ <- NA
                        prozent_fuellung <- NA
                      }
                    }else{
                      dbh.circle <- NA
                      x.circle <- NA
                      y.circle <- NA
                      anteil_punkte_circ <- NA
                      prozent_fuellung <- NA
                    }
                    
                    #Daten holen die nicht vervoxelt sind - nur wenn wir Zeit haben - max 1-3 Kerne
                    # voll_data <- thin3@data
                    
                    voll_data <- filter_poi(sliVox, Z >= u.grenzen.vec[j] & Z<=u.grenzen.vec[j]+z.breite) #das ist die Alternative
                    voll_data <- data.frame("X" = voll_data$X, "Y" = voll_data$Y, "Z" = voll_data$Z,
                                            "cluster" = voll_data$cluster, "Intensity" = voll_data$Intensity)
                    
                    
                    #voll_data <- voll_data[voll_data$Z>=u.grenzen.vec[j]&voll_data$Z<=u.grenzen.vec[j]+z.breite, ]
                    if(!is.na(x.circle)){
                      gross_circ <- point.in.polygon(voll_data$X, voll_data$Y, circle.vert.gross[, 1], circle.vert.gross[, 2])
                      plot.j <- voll_data[which(gross_circ==1), ]
                    }else{
                      if(q<0.75){
                        x.max <- max(plot.j$X)
                        x.min <- min(plot.j$X)
                        y.max <- max(plot.j$Y)
                        y.min <- min(plot.j$Y)
                        plot.j <- voll_data[voll_data$X>=x.min & voll_data$X<=x.max & voll_data$Y>=y.min & voll_data$Y<=y.max, ]
                      }else{
                        #Zuschlag fuer Ausschneiden
                        if(dbh.ellipse<30){zuschlag=0.06}else{zuschlag=0.09}
                        ellipse.vert <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3]+zuschlag, b=geom.ellipse[4]+zuschlag,
                                                         angle=180/pi*geom.ellipse[5], steps=n.vertices)
                        gross_ell <- point.in.polygon(voll_data$X, voll_data$Y, ellipse.vert[, 1], ellipse.vert[, 2])
                        plot.j <- voll_data[which(gross_ell==1), ]
                      }
                    }
                    
                    #anzahl.punkte <- nrow(plot.j)
                    if(allFiles) plot(plot.j$X, plot.j$Y, asp=1)
                    
                    #Endgueltige Ellipse auf nicht vervoxelte Daten
                    par.ellipse <- EllipseDirectFit(cbind(plot.j$X, plot.j$Y))
                    geom.ellipse <- as.vector(AtoG(par.ellipse)$ParG)
                    n.vertices <- 400
                    ellipse.vert <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3], b=geom.ellipse[4],
                                                     angle=180/pi*geom.ellipse[5], steps=n.vertices)
                    q <- geom.ellipse[4]/geom.ellipse[3]
                    if(q>1){q <- geom.ellipse[3]/geom.ellipse[4]}
                    x.ellipse <- geom.ellipse[1]
                    y.ellipse <- geom.ellipse[2]
                    dbh.ellipse <- sqrt(((2*geom.ellipse[4])^2+(2*geom.ellipse[3])^2)/2) * 100
                    if(allFiles) points(ellipse.vert, type="l", lwd=3, lty=2, col=4) #Ellipse zeichnen
                    
                    #Endgueltiger Kreis auf nicht vevoxelte Daten mit anderem Paket
                    BHD_f_fit <- if(is.na(dbh.circle)){
                      x.start <- x.ellipse
                      y.start <- y.ellipse
                      r.start <- dbh.ellipse/100/2
                    }else{
                      x.start <- x.circle
                      y.start <- y.circle
                      r.start <- dbh.circle/100/2
                    }
                    
                    dbh.circle.2 <- r.start*100*2
                    tryCatch({
                      BHD_fit <- LMcircleFit(plot.j[, 1:2], ParIni=c(x.start, y.start, r.start),  IterMAX = 20)
                      pos.dbh2.x <- BHD_fit[1]
                      pos.dbh2.y <- BHD_fit[2]
                      dbh.circle.2 <- BHD_fit[3]*2*100
                      if(allFiles) circle(pos.dbh2.x, pos.dbh2.y, r=dbh.circle.2/100/2, border=3, lwd=2) #neuer Kreis-Fit
                    }, error = function(error_condition){
                    })
                    
                    # dbh.circle.pratt <- CircleFitByPratt(plot.j[, 1:2])[3]*2*100 #Anderer Circ algo
                    
                    if(allFiles){
                      if(!is.na(dbh.circle)){
                        circle(x.circle, y.circle, r=dbh.circle/100/2, border=2, lwd=2) #alter Kreis-Fit
                      }
                      legend("topleft", legend=c("alter Circ", "neuer Circ", "neue Ell"), lty=c(1, 1, 2), col=c(2, 3, 4), lwd=2)
                    }
                    
                    flaeche <- NA
                    
                    #polygon(out2$x[pathX, ], border="red")
                    
                    # zentrum <- as.data.frame(centroid(out2$x[pathX, ]))
                    # colnames(zentrum) <- c("x", "y")
                    #points(zentrum$x, zentrum$y, col=3, pch=15, cex=3)
                    
                    #Ellipse einzeichnen
                    
                    #plot(plot.j$X, plot.j$Y, asp=1) #Zeichnen der Punkte
                    schoener.kreis <- NA
                    
                    if(allFiles) title(main = paste(clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ", ", u.grenzen.vec[j], ", flaeche:", round(flaeche, 3),
                                                    ", \n d.ell:", round(dbh.ellipse, 2), " q.ell:", round(q, 3), " d.circ2:", round(dbh.circle.2, 2), " d.circ:", round(dbh.circle, 2), " n.circ:",
                                                    round(anteil_punkte_circ, 3), " n.punkt:", anzahl.punkte, " prozFuell:", round(prozent_fuellung, 1), sep=""))
                    
                    if(allFiles) dev.off()
                    
                    if(allFiles){
                      png(paste(path.output.cluster.endgraph, "cluster_", clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], "_", "hoehe_", u.grenzen.vec[j], "_gam.png", sep=""),
                          type = "cairo", width = 1200, height = 1200)
                      
                      
                      par(mfrow=c(2, 1))
                      
                    }
                    #plot(plot.j$X, plot.j$Y, asp=1)
                    auffalten <- plot.j
                    zentrum.x <- ifelse(is.na(x.circle), x.ellipse, x.circle)
                    zentrum.y <- ifelse(is.na(y.circle), y.ellipse, y.circle)
                    
                    x_null <- auffalten$X - zentrum.x
                    y_null <- auffalten$Y - zentrum.y
                    
                    auffalten$dist <- sqrt(x_null**2+y_null**2)
                    auffalten$winkel <- atan2(y_null, x_null)
                    if(allFiles) plot(auffalten$winkel, auffalten$dist, xlim=c(-pi, pi), ylim=c(0, 0.5), xlab="Winkel", ylab="Distanz")
                    
                    
                    # plot(auffalten$winkel, auffalten$dist, xlim=c(0, 360), ylim=c(0, 0.5), xlab="Winkel", ylab="Distanz")
                    gam_winkel <- gam(dist ~ s(winkel, bs="cc"), data = auffalten)
                    #gam.check(gam_winkel)
                    summary(gam_winkel)
                    res_gam <- residuals(gam_winkel) #Residuen vom GAM
                    #hist(res_gam)
                    sd(res_gam) #Standardabweichung Residuen #0.02072122
                    
                    st_res_gam <- (res_gam - mean(res_gam))/sd(res_gam) #Standardisierte Residuen vom Gam
                    #hist(st_res_gam)
                    
                    
                    gamfunc <- function(winkel) { predict(gam_winkel, newdata=data.frame("winkel"=winkel)) }
                    
                    werte <- data.frame("winkel"=seq(-pi, pi, 2*pi/360))
                    werte$ergeb <- gamfunc(werte$winkel)
                    
                    if(allFiles) lines(werte$winkel, werte$ergeb, col=2, lwd=2)
                    
                    if(allFiles) title(main = paste(clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ", ", u.grenzen.vec[j], sep=""))
                    
                    #?gam
                    
                    vorhersage <- data.frame("winkel"=seq(-pi, pi, 2*pi/360))
                    vorhersage$dist <- gamfunc(vorhersage$winkel)
                    vorhersage$X <- zentrum.x + cos(vorhersage$winkel) * vorhersage$dist
                    vorhersage$Y <- zentrum.y + sin(vorhersage$winkel) * vorhersage$dist
                    
                    if(allFiles) plot(plot.j$X, plot.j$Y, asp=1) #Zeichnen der Punkte
                    if(allFiles) points(vorhersage$X, vorhersage$Y, col=2, cex=1.5, pch=18)
                    
                    w <- owin(poly=list(x=c(vorhersage$X), y=c(vorhersage$Y)))
                    if(allFiles) plot(w, add=T)
                    
                    #d aus flaeche
                    area <- area.owin(w)
                    d.gam <- sqrt(area/pi*4) * 100 #d aus flaeche
                    
                    if(allFiles) title(main = paste(clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ", ", u.grenzen.vec[j], ", \n d.circ: ",
                                                    round(dbh.circle, 2), ", d.circ2: ", round(dbh.circle.2, 2), ", d.ell: ", round(dbh.ellipse, 2), ", d.gam: ", round(d.gam, 2), sep=""))
                    
                    d.gam
                    dbh.circle
                    dbh.circle.2
                    dbh.ellipse
                    
                    if(allFiles) dev.off()
                    
                    # axis(side=1, at=r.rel, labels = round(r.rel, 3), col.axis="red", col.ticks = "red")
                    
                    #kreis <- data.frame("cluster"=clusterIndex, "x"=x.circle, "y"=y.circle, "BHD" =dbh.circle.1)
                    kreis <- data.frame("cluster"=clusterIndex, "cluster2"=cluster2[clust2], "cluster3"=cluster3[clust3], "z.unten"=u.grenzen.vec[j],
                                        "flaeche"=flaeche, "x.ell" = x.ellipse, "y.ell" = y.ellipse, "d.ell" = dbh.ellipse,
                                        "q.ell" = q, "sk"= schoener.kreis,  "x.circ"= x.circle, "y.circ"=y.circle, "d.circ" = dbh.circle,
                                        "d.circ2" = dbh.circle.2, "n.circ" = anteil_punkte_circ, "n.punkt"= anzahl.punkte,
                                        "prozFuell" = prozent_fuellung, "d.gam"=d.gam)
                    
                    #Qgam - hier kein qgam
                    
                  }else{
                    kreis <- data.frame("cluster"=clusterIndex, "cluster2"=cluster2[clust2], "cluster3"=cluster3[clust3], "z.unten"=u.grenzen.vec[j],
                                        "flaeche"=NA, "x.ell" = NA, "y.ell" = NA, "d.ell" = NA,
                                        "q.ell" = NA, "sk"= NA,  "x.circ"= NA, "y.circ"=NA, "d.circ" = NA, "d.circ2" = NA,
                                        "n.circ" = NA, "n.punkt"= NA, "prozFuell" = NA, "d.gam" = NA)
                    
                  }
                  
                }else{
                  kreis <- data.frame("cluster"=clusterIndex, "cluster2"=cluster2[clust2], "cluster3"=cluster3[clust3], "z.unten"=u.grenzen.vec[j],
                                      "flaeche"=NA, "x.ell" = NA, "y.ell" = NA, "d.ell" = NA,
                                      "q.ell" = NA, "sk"= NA,  "x.circ"= NA, "y.circ"=NA, "d.circ" = NA, "d.circ2" = NA,
                                      "n.circ" = NA, "n.punkt"= NA, "prozFuell" = NA, "d.gam" = NA)
                  
                }
                
                
              }else{ #if(nrow(plot.i)>1300){
                anzahl.punkte <- nrow(plot.i)
                
                plot.j <- plot.i
                # ?circMclust
                if(allFiles){
                  png(paste(path.output.cluster.endgraph, "cluster_", clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], "_", "hoehe_", u.grenzen.vec[j] , ".png", sep=""),
                      type = "cairo", width = 1000, height = 1000)
                  par(mfrow=c(1, 1))
                }
                #plot(plot.j$X, plot.j$Y, asp=1) #Zeichnen der Punkte
                
                #Endgueltige Ellipse
                par.ellipse <- EllipseDirectFit(cbind(plot.j$X, plot.j$Y))
                geom.ellipse <- as.vector(AtoG(par.ellipse)$ParG)
                n.vertices <- 400
                ellipse.vert <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3], b=geom.ellipse[4],
                                                 angle=180/pi*geom.ellipse[5], steps=n.vertices)
                q <- geom.ellipse[4]/geom.ellipse[3]
                if(q>1){q <- geom.ellipse[3]/geom.ellipse[4]}
                x.ellipse <- geom.ellipse[1]
                y.ellipse <- geom.ellipse[2]
                dbh.ellipse <- sqrt(((2*geom.ellipse[4])^2+(2*geom.ellipse[3])^2)/2) * 100
                
                
                #Endgueltiger Kreis
                #Circle Fit
                set.seed(15)
                if( (class(try(
                  if(fast){
                    find.clust.circ <- circMclust( datax=plot.j$X,
                                                   datay=plot.j$Y,
                                                   method="const",
                                                   # nx=1, ny=1, nr=1, # nx=10, ny=10, nr=25,
                                                   nx=8, ny=8, nr=4, #nx=25, ny=25, nr=5,
                                                   minsr=0.01, maxsr=0.5,
                                                   nc=1, # nc=1,
                                                   # minsd=0.1, maxsd=0.5,
                                                   bw=0.01 # bw=0.05
                                                   
                    )
                  } else {
                    # extended nx, ny and nr, takes way longer! would need more extensive testing, but before other efficiency measures
                    find.clust.circ <- circMclust( datax=plot.j$X,
                                                   datay=plot.j$Y,
                                                   method="const",
                                                   # nx=1, ny=1, nr=1, # nx=10, ny=10, nr=25,
                                                   nx=10, ny=10, nr=5, # nx=25, ny=25, nr=5, MAX DETAIL
                                                   minsr=0.01, maxsr=0.5,
                                                   nc=1, # nc=1,
                                                   # minsd=0.1, maxsd=0.5,
                                                   bw=0.01 # bw=0.05
                    )
                  }
                ))=="try-error")==FALSE ){
                  try.error <- FALSE
                }else{
                  try.error <- TRUE
                }
                
                
                if(is.null(find.clust.circ)&try.error==FALSE){ #Problem beim Circle Fit abfangen
                  dbh.circle <- NA
                  x.circle <- NA
                  y.circle <- NA
                  anteil_punkte_circ <- NA
                }
                
                if(!is.null(find.clust.circ)&try.error==F){ #Problem beim Circle Fit abfangen
                  if(colnames(find.clust.circ)[4]=="value"){ #Problem beim Circle Fit abfangen
                    if(is.na(as.numeric(bestMclust(find.clust.circ, 1))[1])==FALSE){ #Problem beim Circle Fit abfangen
                      par.circle.1 <- as.numeric(bestMclust(find.clust.circ, 1))
                      dbh.circle <- round(2 * par.circle.1[3] * 100, 6)
                      x.circle <- par.circle.1[1]
                      y.circle <- par.circle.1[2]
                      
                      #points(par.circle.1[1], par.circle.1[2], cex=1.5, col=6, pch=13)
                      #circle(par.circle.1[1], par.circle.1[2], par.circle.1[3], border=6, lwd=2)
                      
                      #Zuschlag fuer Ausschneiden
                      if(dbh.circle<30){zuschlag=0.06}else{zuschlag=0.09}
                      
                      circle.vert.gross <- calculateCircle(par.circle.1[1], par.circle.1[2], par.circle.1[3]+zuschlag, steps=n.vertices)
                      circle.vert.klein <- calculateCircle(par.circle.1[1], par.circle.1[2], par.circle.1[3]*0.8, steps=n.vertices)
                      
                      klein_circ <- point.in.polygon(plot.j$X, plot.j$Y, circle.vert.klein[, 1], circle.vert.klein[, 2])
                      gross_circ <- point.in.polygon(plot.j$X, plot.j$Y, circle.vert.gross[, 1], circle.vert.gross[, 2])
                      
                      # points(circle.vert.klein, type="l", lwd=3, lty=2, col=2)
                      # points(circle.vert.gross, type="l", lwd=3, lty=2, col=3)
                      
                      anteil_punkte_circ <- sum(klein_circ)/sum(gross_circ)
                      
                      kreis.segment <- plot.j[gross_circ==1&klein_circ!=1, ]
                      #points(kreis.segment$X, kreis.segment$Y, asp=1, col=2)
                      kreis.segment$winkel <- 1
                      neuew=1
                      for(neuew in 1:nrow(kreis.segment)){
                        kreis.segment[neuew, "winkel"] <- angle_points(kreis.segment[neuew, "X"], kreis.segment[neuew, "Y"],
                                                                       x.circle, y.circle)
                      }
                      
                      #hist(kreis.segment$winkel, breaks=seq(0, 350, 10))
                      winkel <- seq(0, 350, 10)
                      breite <- 10
                      wink=1
                      
                      for(wink in 1:length(winkel)){
                        untere_grenz <- winkel[wink]
                        obere_grenz <- winkel[wink]+breite
                        sub.i <- kreis.segment[kreis.segment$winkel>=untere_grenz&kreis.segment$winkel<obere_grenz, ]
                        if(!nrow(sub.i)==0){
                          gefuellt <- 1
                        }else{
                          gefuellt <- 0
                        }
                        if(wink==1){
                          gefuellt.out <- gefuellt
                        }else{
                          gefuellt.out <- c(gefuellt.out, gefuellt)}
                      }
                      
                      prozent_fuellung <- sum(gefuellt.out)/length(gefuellt.out)*100
                      
                      
                    }else{
                      dbh.circle <- NA
                      x.circle <- NA
                      y.circle <- NA
                      anteil_punkte_circ <- NA
                      prozent_fuellung <- NA
                    }
                    
                  }else{
                    dbh.circle <- NA
                    x.circle <- NA
                    y.circle <- NA
                    anteil_punkte_circ <- NA
                    prozent_fuellung <- NA
                  }
                }else{
                  dbh.circle <- NA
                  x.circle <- NA
                  y.circle <- NA
                  anteil_punkte_circ <- NA
                  prozent_fuellung <- NA
                }
                
                #Daten holen die nicht vervoxelt sind - nur wenn wir Zeit haben - max 1-3 Kerne
                # voll_data <- thin3@data #nicht vervoxelt
                
                #voll_data <- tab.neu #das ist die Alternative - vervoxelt
                
                #voll_data <- voll_data[voll_data$Z>=u.grenzen.vec[j]&voll_data$Z<=u.grenzen.vec[j]+z.breite, ]
                voll_data <- filter_poi(sliVox, Z >= u.grenzen.vec[j] & Z<=u.grenzen.vec[j]+z.breite) #das ist die Alternative
                voll_data <- data.frame("X" = voll_data$X, "Y" = voll_data$Y, "Z" = voll_data$Z,
                                        "cluster" = voll_data$cluster, "Intensity" = voll_data$Intensity)
                
                if(!is.na(x.circle)){
                  gross_circ <- point.in.polygon(voll_data$X, voll_data$Y, circle.vert.gross[, 1], circle.vert.gross[, 2])
                  plot.j <- voll_data[which(gross_circ==1), ]
                }else{
                  if(q<0.75){
                    x.max <- max(plot.j$X)
                    x.min <- min(plot.j$X)
                    y.max <- max(plot.j$Y)
                    y.min <- min(plot.j$Y)
                    plot.j <- voll_data[voll_data$X>=x.min & voll_data$X<=x.max & voll_data$Y>=y.min & voll_data$Y<=y.max, ]
                  }else{
                    #Zuschlag fuer Ausschneiden
                    if(dbh.ellipse<30){zuschlag=0.06}else{zuschlag=0.09}
                    ellipse.vert <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3]+zuschlag, b=geom.ellipse[4]+zuschlag,
                                                     angle=180/pi*geom.ellipse[5], steps=n.vertices)
                    gross_ell <- point.in.polygon(voll_data$X, voll_data$Y, ellipse.vert[, 1], ellipse.vert[, 2])
                    plot.j <- voll_data[which(gross_ell==1), ]
                  }
                }
                
                #anzahl.punkte <- nrow(plot.j)
                if(allFiles) plot(plot.j$X, plot.j$Y, asp=1)
                
                #Endgueltige Ellipse auf nicht vervoxelte Daten
                par.ellipse <- EllipseDirectFit(cbind(plot.j$X, plot.j$Y))
                geom.ellipse <- as.vector(AtoG(par.ellipse)$ParG)
                n.vertices <- 400
                ellipse.vert <- calculateEllipse(x=geom.ellipse[1], y=geom.ellipse[2], a=geom.ellipse[3], b=geom.ellipse[4],
                                                 angle=180/pi*geom.ellipse[5], steps=n.vertices)
                q <- geom.ellipse[4]/geom.ellipse[3]
                if(q>1){q <- geom.ellipse[3]/geom.ellipse[4]}
                x.ellipse <- geom.ellipse[1]
                y.ellipse <- geom.ellipse[2]
                dbh.ellipse <- sqrt(((2*geom.ellipse[4])^2+(2*geom.ellipse[3])^2)/2) * 100
                if(allFiles) points(ellipse.vert, type="l", lwd=3, lty=2, col=4) #Ellipse zeichnen
                
                #Endgueltiger Kreis auf nicht vevoxelte Daten mit anderem Paket
                BHD_f_fit <- if(is.na(dbh.circle)){
                  x.start <- x.ellipse
                  y.start <- y.ellipse
                  r.start <- dbh.ellipse/100/2
                }else{
                  x.start <- x.circle
                  y.start <- y.circle
                  r.start <- dbh.circle/100/2
                }
                
                
                dbh.circle.2 <- r.start*100*2
                tryCatch({
                  BHD_fit <- LMcircleFit(plot.j[, 1:2], ParIni=c(x.start, y.start, r.start),  IterMAX = 20)
                  pos.dbh2.x <- BHD_fit[1]
                  pos.dbh2.y <- BHD_fit[2]
                  dbh.circle.2 <- BHD_fit[3]*2*100
                  if(allFiles) circle(pos.dbh2.x, pos.dbh2.y, r=dbh.circle.2/100/2, border=3, lwd=2) #neuer Kreis-Fit
                }, error = function(error_condition){
                })
                
                #dbh.circle.pratt <- CircleFitByPratt(plot.j[, 1:2])[3]*2*100 #Anderer Circ algo
                if(allFiles) {
                  if(!is.na(dbh.circle)){
                    circle(x.circle, y.circle, r=dbh.circle/100/2, border=2, lwd=2) #alter Kreis-Fit
                  }
                  legend("topleft", legend=c("alter Circ", "neuer Circ", "neue Ell"), lty=c(1, 1, 2), col=c(2, 3, 4), lwd=2)
                }
                
                flaeche <- NA
                #polygon(out2$x[pathX, ], border="red")
                
                # zentrum <- as.data.frame(centroid(out2$x[pathX, ]))
                # colnames(zentrum) <- c("x", "y")
                #points(zentrum$x, zentrum$y, col=3, pch=15, cex=3)
                
                #Ellipse einzeichnen
                
                #plot(plot.j$X, plot.j$Y, asp=1) #Zeichnen der Punkte
                schoener.kreis <- NA
                
                if(allFiles) {
                  title(main = paste(clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ", ", u.grenzen.vec[j], ", flaeche:", round(flaeche, 3),
                                     ", \n d.ell:", round(dbh.ellipse, 2), " q.ell:", round(q, 3), " d.circ2:", round(dbh.circle.2, 2), " d.circ:", round(dbh.circle, 2), " n.circ:",
                                     round(anteil_punkte_circ, 3), " n.punkt:", anzahl.punkte, " prozFuell:", round(prozent_fuellung, 1), sep=""))
                  
                  dev.off()
                  
                  png(paste(path.output.cluster.endgraph, "cluster_", clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], "_", "hoehe_", u.grenzen.vec[j], "_gam.png", sep=""),
                      type = "cairo", width = 1200, height = 1200)
                  
                  par(mfrow=c(2, 1))
                }
                
                
                
                
                #plot(plot.j$X, plot.j$Y, asp=1)
                auffalten <- plot.j
                zentrum.x <- ifelse(is.na(x.circle), x.ellipse, x.circle)
                zentrum.y <- ifelse(is.na(y.circle), y.ellipse, y.circle)
                
                x_null <- auffalten$X - zentrum.x
                y_null <- auffalten$Y - zentrum.y
                
                auffalten$dist <- sqrt(x_null**2+y_null**2)
                auffalten$winkel <- atan2(y_null, x_null)
                if(allFiles) plot(auffalten$winkel, auffalten$dist, xlim=c(-pi, pi), ylim=c(0, 0.5), xlab="Winkel", ylab="Distanz")
                
                # plot(auffalten$winkel, auffalten$dist, xlim=c(0, 360), ylim=c(0, 0.5), xlab="Winkel", ylab="Distanz")
                gam_winkel <- gam(dist ~ s(winkel, bs="cc"), data = auffalten)
                #gam.check(gam_winkel)
                summary(gam_winkel)
                res_gam <- residuals(gam_winkel) #Residuen vom GAM
                #hist(res_gam)
                sd(res_gam) #Standardabweichung Residuen #0.02072122
                
                st_res_gam <- (res_gam - mean(res_gam))/sd(res_gam) #Standardisierte Residuen vom Gam
                #hist(st_res_gam)
                
                
                gamfunc <- function(winkel) { predict(gam_winkel, newdata=data.frame("winkel"=winkel)) }
                
                werte <- data.frame("winkel"=seq(-pi, pi, 2*pi/360))
                werte$ergeb <- gamfunc(werte$winkel)
                
                if(allFiles){
                  lines(werte$winkel, werte$ergeb, col=2, lwd=2)
                  
                  title(main = paste(clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ", ", u.grenzen.vec[j], sep=""))
                }
                #?gam
                
                vorhersage <- data.frame("winkel"=seq(-pi, pi, 2*pi/360))
                vorhersage$dist <- gamfunc(vorhersage$winkel)
                vorhersage$X <- zentrum.x + cos(vorhersage$winkel) * vorhersage$dist
                vorhersage$Y <- zentrum.y + sin(vorhersage$winkel) * vorhersage$dist
                
                
                if(allFiles){
                  plot(plot.j$X, plot.j$Y, asp=1) #Zeichnen der Punkte
                  points(vorhersage$X, vorhersage$Y, col=2, cex=1.5, pch=18)
                }
                
                #w <- owin(poly=list(x=c(rev(vorhersage$X)), y=c(rev(vorhersage$Y))))
                w <- owin(poly=list(x=c(vorhersage$X), y=c(vorhersage$Y)))
                if(allFiles) plot(w, add=T)
                
                
                #d aus flaeche
                area <- area.owin(w)
                d.gam <- sqrt(area/pi*4) * 100 #d aus flaeche
                
                if(allFiles){
                  title(main = paste(clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ", ", u.grenzen.vec[j], ", \n d.circ: ",
                                     round(dbh.circle, 2), ", d.circ2: ", round(dbh.circle.2, 2), ", d.ell: ", round(dbh.ellipse, 2), ", d.gam: ", round(d.gam, 2), sep=""))
                  
                }
                d.gam
                dbh.circle
                dbh.circle.2
                dbh.ellipse
                
                if(allFiles){
                  dev.off()
                }
                
                # axis(side=1, at=r.rel, labels = round(r.rel, 3), col.axis="red", col.ticks = "red")
                
                #kreis <- data.frame("cluster"=clusterIndex, "x"=x.circle, "y"=y.circle, "BHD" =dbh.circle.1)
                kreis <- data.frame("cluster"=clusterIndex, "cluster2"=cluster2[clust2], "cluster3"=cluster3[clust3], "z.unten"=u.grenzen.vec[j],
                                    "flaeche"=flaeche, "x.ell" = x.ellipse, "y.ell" = y.ellipse, "d.ell" = dbh.ellipse,
                                    "q.ell" = q, "sk"= schoener.kreis,  "x.circ"= x.circle, "y.circ"=y.circle, "d.circ" = dbh.circle,
                                    "d.circ2" = dbh.circle.2, "n.circ" = anteil_punkte_circ, "n.punkt"= anzahl.punkte,
                                    "prozFuell" = prozent_fuellung, "d.gam"=d.gam)
                
                #Qgam - hier kein qgam
              }
              
            }else{
              kreis <- data.frame("cluster"=clusterIndex, "cluster2"=cluster2[clust2], "cluster3"=cluster3[clust3], "z.unten"=u.grenzen.vec[j],
                                  "flaeche"=NA, "x.ell" = NA, "y.ell" = NA, "d.ell" = NA,
                                  "q.ell" = NA, "sk"= NA,  "x.circ"= NA, "y.circ"=NA, "d.circ" = NA, "d.circ2" = NA,
                                  "n.circ" = NA, "n.punkt"= NA, "prozFuell" = NA, "d.gam" = NA)
              
            }
            
          }else{
            kreis <- data.frame("cluster"=clusterIndex, "cluster2"=cluster2[clust2], "cluster3"=cluster3[clust3], "z.unten"=u.grenzen.vec[j],
                                "flaeche"=NA, "x.ell" = NA, "y.ell" = NA, "d.ell" = NA,
                                "q.ell" = NA, "sk"= NA,  "x.circ"= NA, "y.circ"=NA, "d.circ" = NA, "d.circ2" = NA,
                                "n.circ" = NA, "n.punkt"= NA, "prozFuell" = NA, "d.gam" = NA)
            
          }
          
          if(is.na(kreis$d.gam)){
            out.3d.data.j <- kreis[, c("cluster", "cluster2", "cluster3")]
            out.3d.data.j$z.schicht.mittel <- u.grenzen.vec[j]+0.125/2
            out.3d.data.j$winkel <- NA
            out.3d.data.j$dist <- NA
            out.3d.data.j$X <- NA
            out.3d.data.j$Y <- NA
            out.3d.data.j$Z <- NA
            out.3d.data.j$x.zentrum <- NA
            out.3d.data.j$y.zentrum <- NA
            out.res.gam.j <- data.frame("cluster"=clusterIndex, "cluster2"=cluster2[clust2], "cluster3"=cluster3[clust3], "z.schicht.mittel" = u.grenzen.vec[j]+0.125/2,
                                        "res_gam"= NA, "st_res_gam" = NA, "inten" = NA)
          } else {
            out.3d.data.j <- auffalten[auffalten$Z>=u.grenzen.vec[j]&auffalten$Z<u.grenzen.vec[j]+0.125, c("X", "Y", "Z", "winkel", "dist")]
            out.3d.data.j$cluster <- clusterIndex
            out.3d.data.j$cluster2 <- cluster2[clust2]
            out.3d.data.j$cluster3 <- cluster3[clust3]
            out.3d.data.j$z.schicht.mittel <- u.grenzen.vec[j]+0.125/2
            out.3d.data.j$x.zentrum <- zentrum.x
            out.3d.data.j$y.zentrum <- zentrum.y
            out.res.gam.j <- data.frame("cluster"=clusterIndex, "cluster2"=cluster2[clust2], "cluster3"=cluster3[clust3], "z.schicht.mittel" = u.grenzen.vec[j]+0.125/2,
                                        "res_gam"= res_gam, "st_res_gam" = st_res_gam, "inten" = plot.j$Intensity)
          }
          
          
          if(j==1){
            out.kreis <- kreis
            out.3d.data <- out.3d.data.j
            out.res <- out.res.gam.j
          }else{
            out.kreis <- rbind(out.kreis, kreis)
            out.3d.data <- rbind(out.3d.data, out.3d.data.j)
            out.res <- rbind(out.res, out.res.gam.j)
          }
        } #for(j in 1:length(u.grenzen.vec)){
        # write.csv2(out.kreis, paste(path.output.cluster.end, "cluster_", clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ".csv", sep=""))
        # write.csv(out.3d.data, paste("F:/Testen_Intensitaet/Ergebnis/tensor_gam/", "cluster_", clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ".csv", sep=""))
        #write.csv2(out.3d.data, paste("F:/Testen_Intensitaet/Ergebnis/", "cluster_", clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ".csv", sep=""))
        
        if(allFiles) write.csv2(out.res,
                                paste(path.output.cluster.residuals, "cluster_",
                                      clusterIndex, "_", cluster2[clust2],
                                      "_", cluster3[clust3], ".csv", sep=""),
                                row.names = F)
        
        
        #Tensorproduktsmooth
        if(nrow(out.3d.data[is.na(out.3d.data$X), ])<9){
          gam_winkel.i <- gam(dist ~ te(winkel, Z, bs=c("cc", "tp")), data = out.3d.data)
          
          # data.pred <- data.frame(unique(cbind(out.3d.data$x.zentrum, out.3d.data$y.zentrum, out.3d.data$z.schicht.mittel)))
          # colnames(data.pred) <- c("zentrum.x", "zentrum.y", "Z")
          #
          # temp <- data.frame(winkel=1:360)
          #
          # data.pred <- merge(data.pred, temp)
          #
          # # plot3d(data.pred$zentrum.x, data.pred$zentrum.y, data.pred$Z, asp=F)
          #
          # data.pred$X.abs <- data.pred$X.rel + data.pred$zentrum.x
          # data.pred$Y.abs <- data.pred$Y.rel + data.pred$zentrum.y
          #
          # {
          #  plot3d(out.3d.data$X, out.3d.data$Y, out.3d.data$Z, asp=F)
          #  points3d(data.pred$X.abs, data.pred$Y.abs, data.pred$Z, col=2, pch=0.5)
          # }
          
          data.pred <- data.frame(winkel=seq(-pi, pi, 2*pi/360), Z=1.3)
          
          data.pred$dist.pred <- predict(gam_winkel.i, newdata=data.pred)
          
          data.pred$X.rel <- cos(data.pred$winkel) * data.pred$dist.pred
          data.pred$Y.rel <- sin(data.pred$winkel) * data.pred$dist.pred
          
          #plot(data.pred$X.rel, data.pred$Y.rel, asp=1)
          
          w <- owin(poly=list(x=c(data.pred$X.rel), y=c(data.pred$Y.rel)))
          #plot(w, add=T)
          
          #d aus flaeche
          area <- area.owin(w)
          d.gam <- sqrt(area/pi*4) * 100 #d aus flaeche
          
          out.kreis$tegam.d_gam <- d.gam
          
          write.csv2(out.kreis, paste(path.output.cluster.end, "cluster_", clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ".csv", sep=""))
        }else{
          
          out.kreis$tegam.d_gam <- NA
          
          write.csv2(out.kreis, paste(path.output.cluster.end, "cluster_", clusterIndex, "_", cluster2[clust2], "_", cluster3[clust3], ".csv", sep=""))
        }
        
        
        #
      } #for(clust3 in 1:lenght(cluster3)){
      
    } #for(clust2 in 1:length(cluster2)){
    
    
  } #if(nrow(plot.clust2)!=0){
  #if(nrow(plot.clust2)!=0){
  
  
}

diameterBeast <- function(fileFinder, dbhPath, ipad = FALSE, allFiles = FALSE, nr_cores = 0,
                          referenced = FALSE, bushPreparation = FALSE, filterSOR = FALSE, fast = TRUE,
                          cutWindow = c(-1000, -1000, 2000), silent = FALSE, dirPath = paste0(getwd(), "/")){
  start <- Sys.time()
  cat("\nStarting diameterBeast()\n")
  cat("Current time is", format(Sys.time(), "%H:%M:%S"), "\n")
  
  
  if(!exists("sliVox")){
    
    cat("Reading in lost slice cluster (normalized):\n")
    # NORMALIZATION
    sliVox <- readLAS(paste0(dbhPath, "slice_cluster_slope.laz"))
    
    tryCatch(
      {
        dtm_z <- raster(paste0(dirPath, v.env$groundPath, fileFinder, "_ground_min.grd"))
        sliVox <- normalize_height(sliVox, dtm_z, na.rm = TRUE)
        cat("done!\n")
      }, error = function(error_condition) {
        dtm_a <- raster(paste0(dirPath, groundPath, fileFinder, "_ground_rough.grd"))
        sliVox <- normalize_height(sliVox, dtm_a, na.rm = TRUE)
      })
    cat("Successfully retrieved normalized slice cluster file!\n\n")
  }
  sliVox <- sliVox
  
  # if you had lost sliVox, read in from paste0(dbhPath, "slice_cluster.laz")
  #sliVoxSafe <- sliVox
  
  XL <- cutWindow[1]
  YL <- cutWindow[2]
  width <- cutWindow[3]
  if(sum(cutWindow == c(-1000, -1000, 2000)) == 3){
    sliVox <- sliVox
    cat("No cutting because of no cutWindow setting!\n")
  } else {
    sliVox <- filter_poi(sliVox, X > XL, X < XL + width, Y > YL, Y < YL + width)
    # CUTTING IS ALREADY DONE ABOVE
  }
  
  
  #tab.neu <- cluster_tab
  cat(fileFinder, " - Processing", length(unique(sliVox@data$cluster)), "clusters",
      "<: from", min(sliVox@data$cluster), "to",  max(sliVox@data$cluster), ":>\n")
  
  
  
  i=1
  j=1
  
  
  
  timePar1 <- Sys.time()
  
  if(nr_cores == 0){
    cat("Automatic number of core setting:\n   ")
    nr_cores <- detectCores()
    cat(nr_cores, "cores available - max ")
    nr_cores <- nr_cores - 1
    
    maxCores <- 10
    cat(maxCores, "cores recommended.\n")
    if(maxCores < 1){
      maxCores <- 1
    }
    if(maxCores <=  nr_cores){
      nr_cores <- maxCores
    }
    #cat("\n")
    # if(las@header@PHB$`Number of point records` > 90000000){
    #   cat("LAS File is too big, reducing back to 5 cores!")
    #   if(nr_cores > 5){
    #     nr_cores <- 5
    #   }
    # }
  }
  
  cluster.vec <- unique(sliVox@data$cluster)
  
  if(nr_cores > length(cluster.vec)){
    nr_cores <- length(cluster.vec)
  }
  
  if(nr_cores == 1){
    
    cat("Going for serial diameter measurement:")
    for(i in 1:length(cluster.vec)){
      if(i%%20 == 1) cat("\n    ")
      cat(paste0("", cluster.vec[i], "-"))
      diameterBeast_i(clusterIndex = cluster.vec[i], dbhPath = dbhPath, fast = fast, allFiles = allFiles)
    }
  }
  else 
  {
    # PARALLEL PART
    cat(fileFinder, " - Going parrrallel on", nr_cores, "cores.\n\n")
    
    library(doParallel)
    cl <- makeCluster(nr_cores)#, outfile="")
    registerDoParallel(cl)
    
    file_parallelProtocol <- paste0(dbhPath, "temp_par_diameterBeast.txt")
    file.create(file_parallelProtocol)
    fdc <<- foreach(i=1:length(cluster.vec),  .errorhandling = 'remove', 
                    .export=c('diameterBeast_i', "sliVox", 'v.env'), 
                    .packages = c("treeX"))%dopar% {
                      t1 <- Sys.time()
                      sink(file_parallelProtocol, append = T)
                      numi <- cluster.vec[i]
                      cat(paste0(numi, "-"))
                      
                      diameterBeast_i(clusterIndex = numi, 
                                      dbhPath = dbhPath, fast = fast, allFiles = allFiles)
                      
                      t2 <- Sys.time()
                      dt <- difftime(t2,t1)
                      dt <- paste(round(dt, 1), units(dt))
                      cat(paste0("+", numi, "in", dt, "\n"))
                      sink()
                      
                    } 
    
    
    stopCluster(cl)
    getDoParWorkers()
    
    #unlink(file_parallelProtocol)
    
    timePar2 <- Sys.time()
    gc()
    
    
    cat("Parallel work is done in a ")
    print.difftime(round(timePar2 - timePar1, 1))
  }
  #}
  stop <- Sys.time()
  cat("Diameter beast is done.\n Total ")
  print.difftime(round(stop - start, 1))
  
  
  
}








####################################################################################################
# 3 FINE CLUSTER CLEARING ##########################################################################
####################################################################################################

#' analyzes, if the measured cluster can be a tree
#'
#' @param circleFile name of the file to be analyzed
fineCluster_i <- function(circleFile){
  
  out <- read.csv2(circleFile)
  out <- out[, -1] #remove row number
  out[out$n.punkt<=90&!is.na(out$n.punkt), c("flaeche", "x.ell", "y.ell", "d.ell", "q.ell", "sk",
                                             "x.circ", "y.circ", "d.circ", "n.circ", "n.punkt", "prozFuell", "d.gam", "d.ob", "d.ob2",
                                             "d.min", "d.max", "d.stip", "tegam.d_ob", "tegam.d_ob2", "tegam.d_min", "tegam.d_max", "tegam.d_stip", "tegam.d_gam")] <- NA #out$n.punkt<=150
  #min punkte auf 90 von 100 heruntergestellt
  out[out$q.ell<=0.60&!is.na(out$x.ell), c("flaeche", "x.ell", "y.ell", "d.ell", "q.ell", "sk")] <- NA #out$q.ell<=0.65
  #Grenze fuer abfalchung Ellipse muss groesser sein wie #0.60 sonst NA
  out[out$prozFuell<25&!is.na(out$x.circ), c("flaeche", "x.ell", "y.ell", "d.ell", "q.ell", "sk",
                                             "x.circ", "y.circ", "d.circ", "n.circ", "n.punkt", "prozFuell", "d.gam", "d.ob", "d.ob2",
                                             "d.min", "d.max", "d.stip", "tegam.d_ob", "tegam.d_ob2", "tegam.d_min", "tegam.d_max", "tegam.d_stip", "tegam.d_gam")] <- NA
  #wir brauchen mind. 25% Fuellung
  
  # out
  # plot(out$z.unten, out$d.circ)
  
  #head(out)
  out.vec <- 1:nrow(out)
  m.vec <- 6 #mit 6 fixiert -> es muessen mind. 6 durch die Schwellenwerte kommen
  
  # Preallocate a list to store results
  output_list <- list()
  index <- 1
  
  for (k in seq_along(m.vec)) {
    #cat("\n", k, "- ")
    
    # Generate all combinations for the current m.vec[k]
    liste <- combn(out.vec, m.vec[k], simplify = FALSE)
    
    # Loop through each combination
    for (l in seq_along(liste)) {
      #cat(l, "-")
      
      # Extract the subset of rows corresponding to the current combination
      subset_data <- out[liste[[l]], ]
      
      # Compute the required statistics
      #sd fuer pos und durchmesser der ellipse
      fertig <- data.frame(
        sdx.ell = sd(subset_data$x.ell),
        sdy.ell = sd(subset_data$y.ell),
        sdbhd.ell = sd(subset_data$d.ell),
        sdx.circ = sd(subset_data$x.circ),
        sdy.circ = sd(subset_data$y.circ),
        sdbhd.circ = sd(subset_data$d.circ),
        #mittlere rel. anzahl an punkten aus den 6
        mean.n.circ = mean(subset_data$n.circ),
        comb.circ = l,
        anzahl.pro = m.vec[k]
      )
      
      # Store the result in the preallocated list
      output_list[[index]] <- fertig
      index <- index + 1
    }
  }
  
  # Combine all results into a single data frame
  output2 <- do.call(rbind, output_list)
  
  grenze.xy <- 0.2 #0.2 passt fuer 44, fuer 68 schon zu gross; ist relativ gross gewaehlt
  grenze.BHD <- 1.85 ##muesste fuer 152 2.2 sein
  #wenn keiner uebrig bleibt, dann gibt es keinen
  
  if(length(unique(output2$mean.n.circ))!=1){
    if(min(output2$mean.n.circ, na.rm = T)<0.35&
       min(output2$sdx.circ, na.rm = T)<0.07&min(output2$sdy.circ, na.rm = T)<0.07){grenze.BHD <- 2.0}
    
    if(min(output2$mean.n.circ, na.rm = T)<0.4&
       min(output2$sdx.circ, na.rm = T)<0.02&min(output2$sdy.circ, na.rm = T)<0.02){grenze.BHD <- 4.7}
    
    if(median(output2$sdbhd.circ, na.rm =T) > 2.7){grenze.BHD <- 0.8}
  }
  
  summary(output2)
  # hist(output2$sdbhd.circ)
  # mean(output2$sdbhd.circ, na.rm =T)
  # median(output2$sdbhd.circ, na.rm =T)
  # skewness(output2$sdbhd.circ, na.rm =T)
  
  #auf Grenzen subsetten
  optim2.circ <- output2[output2$sdx.circ<=grenze.xy&output2$sdy.circ<=grenze.xy&output2$sdbhd.circ<=grenze.BHD&!is.na(output2$sdbhd.circ), ]
  optim2.ell <- output2[output2$sdx.ell<=grenze.xy&output2$sdy.ell<=grenze.xy&output2$sdbhd.ell<=grenze.BHD&!is.na(output2$sdbhd.ell), ]
  
  if(nrow(optim2.circ)!=0){ #wenn im Kreis was drinnen ist
    optim3.circ <- optim2.circ #
    
    ende.circ <- optim3.circ[optim3.circ$sdbhd.circ<=quantile(optim3.circ$sdbhd.circ, 0.05), ] #nimm die 5% der genauesten hinsichtlich des bhd
    
    liste6 <- combn(out.vec, 6, simplify = FALSE)
    liste6 <- lapply(liste6, as.character)
    vec.6 <- which(sapply(liste6, FUN=function(X) "3" %in% X)) #weil 3 der auf 1.250 BHD-Hoehe ist #Schauen ob die in den Kombis vorhanden ist
    
    ende2.circ <- ende.circ[ende.circ$comb.circ %in% vec.6, ] #
    
    if(nrow(ende2.circ)!=0){#wenn ein Circle vorhanden ist, nehme den besten, hinsichtlich des sdBHD
      ende3.circ <- ende2.circ[ende2.circ$sdbhd.circ==min(ende2.circ$sdbhd.circ), ][1, ]
      
      liste6 <- combn(out.vec, 6, simplify = FALSE)
      #schauen was gebraucht wird - rest NA machen
      wahl_circ <- out[3, ]
      wahl_circ[, c(6:9)] <- NA
      wahl_circ6 <- out[liste6[[ende3.circ$comb.circ]], ]
      wahl_circ6[, c(6:9)] <- NA
    }else{ende3.circ <- ende2.circ}
    
    if(nrow(ende3.circ)==0){ #wenn es fuer diesen Druchmesser keine Variante gibt, dann anderen suchen
      
      vec.6 <- which(sapply(liste6, FUN=function(X) "2" %in% X)) #2 ist der unter dem BHD auf 1.125
      
      ende2.circ <- ende.circ[ende.circ$comb.circ %in% vec.6, ]
      
      if(nrow(ende2.circ)!=0){
        ende3.circ <- ende2.circ[ende2.circ$sdbhd.circ==min(ende2.circ$sdbhd.circ), ][1, ]
        
        liste6 <- combn(out.vec, 6, simplify = FALSE)
        #schauen was gebraucht wird - rest NA machen
        wahl_circ <- out[2, ]
        wahl_circ[, c(6:9)] <- NA
        wahl_circ6 <- out[liste6[[ende3.circ$comb.circ]], ]
        wahl_circ6[, c(6:9)] <- NA
      }else{ende3.circ <- ende2.circ}
      
      if(nrow(ende3.circ)==0){
        
        vec.6 <- which(sapply(liste6, FUN=function(X) "1" %in% X)) #
        
        ende2.circ <- ende.circ[ende.circ$comb.circ %in% vec.6, ]
        
        if(nrow(ende2.circ)!=0){
          ende3.circ <- ende2.circ[ende2.circ$sdbhd.circ==min(ende2.circ$sdbhd.circ), ][1, ]
          
          liste6 <- combn(out.vec, 6, simplify = FALSE)
          #schauen was gebraucht wird - rest NA machen
          wahl_circ <- out[1, ]
          wahl_circ[, c(6:9)] <- NA
          wahl_circ6 <- out[liste6[[ende3.circ$comb.circ]], ]
          wahl_circ6[, c(6:9)] <- NA
        }else{ende3.circ <- ende2.circ}
        
        if(nrow(ende3.circ)==0){
          
          vec.6 <- which(sapply(liste6, FUN=function(X) "4" %in% X)) #
          
          ende2.circ <- ende.circ[ende.circ$comb.circ %in% vec.6, ]
          
          if(nrow(ende2.circ)!=0){
            ende3.circ <- ende2.circ[ende2.circ$sdbhd.circ==min(ende2.circ$sdbhd.circ), ][1, ]
            
            liste6 <- combn(out.vec, 6, simplify = FALSE)
            #schauen was gebraucht wird - rest NA machen
            wahl_circ <- out[4, ]
            wahl_circ[, c(6:9)] <- NA
            wahl_circ6 <- out[liste6[[ende3.circ$comb.circ]], ]
            wahl_circ6[, c(6:9)] <- NA
          }else{ende3.circ <- ende2.circ}
          
          if(nrow(ende3.circ)==0){
            
            vec.6 <- which(sapply(liste6, FUN=function(X) "5" %in% X)) #
            
            ende2.circ <- ende.circ[ende.circ$comb.circ %in% vec.6, ]
            
            if(nrow(ende2.circ)!=0){
              ende3.circ <- ende2.circ[ende2.circ$sdbhd.circ==min(ende2.circ$sdbhd.circ), ][1, ]
              
              liste6 <- combn(out.vec, 6, simplify = FALSE)
              #schauen was gebraucht wird - rest NA machen
              wahl_circ <- out[5, ]
              wahl_circ[, c(6:9)] <- NA
              wahl_circ6 <- out[liste6[[ende3.circ$comb.circ]], ]
              wahl_circ6[, c(6:9)] <- NA
            }else{ende3.circ <- ende2.circ}
            
            if(nrow(ende3.circ)==0){
              
              vec.6 <- which(sapply(liste6, FUN=function(X) "6" %in% X)) #
              
              ende2.circ <- ende.circ[ende.circ$comb.circ %in% vec.6, ]
              
              if(nrow(ende2.circ)!=0){
                ende3.circ <- ende2.circ[ende2.circ$sdbhd.circ==min(ende2.circ$sdbhd.circ), ][1, ]
                
                liste6 <- combn(out.vec, 6, simplify = FALSE)
                #schauen was gebraucht wird - rest NA machen
                wahl_circ <- out[6, ]
                wahl_circ[, c(6:9)] <- NA
                wahl_circ6 <- out[liste6[[ende3.circ$comb.circ]], ]
                wahl_circ6[, c(6:9)] <- NA
              }else{ende3.circ <- ende2.circ}
              
              if(nrow(ende3.circ)==0){
                
                vec.6 <- which(sapply(liste6, FUN=function(X) "7" %in% X)) #
                
                ende2.circ <- ende.circ[ende.circ$comb.circ %in% vec.6, ]
                
                if(nrow(ende2.circ)!=0){
                  ende3.circ <- ende2.circ[ende2.circ$sdbhd.circ==min(ende2.circ$sdbhd.circ), ][1, ]
                  
                  liste6 <- combn(out.vec, 6, simplify = FALSE)
                  #schauen was gebraucht wird - rest NA machen
                  wahl_circ <- out[7, ]
                  wahl_circ[, c(6:9)] <- NA
                  wahl_circ6 <- out[liste6[[ende3.circ$comb.circ]], ]
                  wahl_circ6[, c(6:9)] <- NA
                }else{ende3.circ <- ende2.circ}
                
                if(nrow(ende3.circ)==0){
                  
                  vec.6 <- which(sapply(liste6, FUN=function(X) "8" %in% X)) #
                  
                  ende2.circ <- ende.circ[ende.circ$comb.circ %in% vec.6, ]
                  
                  if(nrow(ende2.circ)!=0){
                    ende3.circ <- ende2.circ[ende2.circ$sdbhd.circ==min(ende2.circ$sdbhd.circ), ][1, ]
                    
                    liste6 <- combn(out.vec, 6, simplify = FALSE)
                    #schauen was gebraucht wird - rest NA machen
                    wahl_circ <- out[8, ]
                    wahl_circ[, c(6:9)] <- NA
                    wahl_circ6 <- out[liste6[[ende3.circ$comb.circ]], ]
                    wahl_circ6[, c(6:9)] <- NA
                  }else{ende3.circ <- ende2.circ}
                  
                  if(nrow(ende3.circ)==0){
                    
                    vec.6 <- which(sapply(liste6, FUN=function(X) "9" %in% X)) #
                    
                    ende2.circ <- ende.circ[ende.circ$comb.circ %in% vec.6, ]
                    
                    if(nrow(ende2.circ)!=0){
                      ende3.circ <- ende2.circ[ende2.circ$sdbhd.circ==min(ende2.circ$sdbhd.circ), ][1, ]
                      
                      liste6 <- combn(out.vec, 6, simplify = FALSE)
                      #schauen was gebraucht wird - rest NA machen
                      wahl_circ <- out[9, ]
                      wahl_circ[, c(6:9)] <- NA
                      wahl_circ6 <- out[liste6[[ende3.circ$comb.circ]], ]
                      wahl_circ6[, c(6:9)] <- NA
                    }else{ende3.circ <- ende2.circ}
                  }
                  
                }
                
              }
              
            }
            
          }
        }
        
      }
      
    }
  }else{
    wahl_circ <- out[1, ]
    wahl_circ[1, 4:ncol(wahl_circ)] <- NA
    wahl_circ6 <- out[1:6, ]
    wahl_circ6[, 4:ncol(wahl_circ6)] <- NA
  }
  
  if(nrow(optim2.ell)!=0){ #wenn in einem was drinnen ist
    optim3.ell <- optim2.ell
    
    if(nrow(optim3.ell)!=0){ #wenn es dann noch was gibt
      
      ende.ell <- optim3.ell[optim3.ell$sdbhd.ell<=quantile(optim3.ell$sdbhd.ell, 0.05), ] #nimm die 5% der genauesten hinsichtlich des bhd
      
      liste6 <- combn(out.vec, 6, simplify = FALSE)
      liste6 <- lapply(liste6, as.character)
      vec.6 <- which(sapply(liste6, FUN=function(X) "3" %in% X)) #weil 3 der auf 1.250 BHD-Hoehe ist #Schauen ob die in den Kombis vorhanden ist
      
      ende2.ell <- ende.ell[ende.ell$comb.circ %in% vec.6, ]
      
      if(nrow(ende2.ell)!=0){#wenn eine Ellipse vorhanden ist, dasselbe
        ende3.ell <- ende2.ell[ende2.ell$sdbhd.ell==min(ende2.ell$sdbhd.ell), ][1, ]
        
        liste6 <- combn(out.vec, 6, simplify = FALSE)
        #schauen was gebraucht wird - rest NA machen
        wahl_ell <- out[3, ]
        wahl_ell[, c(11:14)] <- NA
        wahl_ell6 <- out[liste6[[ende3.ell$comb.circ]], ]
        wahl_ell6[, c(11:14)] <- NA
      }else{ende3.ell <- ende2.ell}
      
      if(nrow(ende3.ell)==0){ #wenn es fuer diesen Druchmesser keine Variante gibt, dann anderen suchen
        
        vec.6 <- which(sapply(liste6, FUN=function(X) "2" %in% X)) #2 ist der unter dem BHD auf 1.125
        
        ende2.ell <- ende.ell[ende.ell$comb.circ %in% vec.6, ]
        
        if(nrow(ende2.ell)!=0){
          ende3.ell <- ende2.ell[ende2.ell$sdbhd.ell==min(ende2.ell$sdbhd.ell), ][1, ]
          
          liste6 <- combn(out.vec, 6, simplify = FALSE)
          #schauen was gebraucht wird - rest NA machen
          wahl_ell <- out[2, ]
          wahl_ell[, c(11:14)] <- NA
          wahl_ell6 <- out[liste6[[ende3.ell$comb.circ]], ]
          wahl_ell6[, c(11:14)] <- NA
        }else{ende3.ell <- ende2.ell}
        
        if(nrow(ende3.ell)==0){
          
          vec.6 <- which(sapply(liste6, FUN=function(X) "1" %in% X)) #
          
          ende2.ell <- ende.ell[ende.ell$comb.circ %in% vec.6, ]
          
          if(nrow(ende2.ell)!=0){
            ende3.ell <- ende2.ell[ende2.ell$sdbhd.ell==min(ende2.ell$sdbhd.ell), ][1, ]
            
            liste6 <- combn(out.vec, 6, simplify = FALSE)
            #schauen was gebraucht wird - rest NA machen
            wahl_ell <- out[1, ]
            wahl_ell[, c(11:14)] <- NA
            wahl_ell6 <- out[liste6[[ende3.ell$comb.circ]], ]
            wahl_ell6[, c(11:14)] <- NA
          }else{ende3.ell <- ende2.ell}
          
          if(nrow(ende3.ell)==0){
            
            vec.6 <- which(sapply(liste6, FUN=function(X) "4" %in% X)) #
            
            ende2.ell <- ende.ell[ende.ell$comb.circ %in% vec.6, ]
            
            if(nrow(ende2.ell)!=0){
              ende3.ell <- ende2.ell[ende2.ell$sdbhd.ell==min(ende2.ell$sdbhd.ell), ][1, ]
              
              liste6 <- combn(out.vec, 6, simplify = FALSE)
              #schauen was gebraucht wird - rest NA machen
              wahl_ell <- out[4, ]
              wahl_ell[, c(11:14)] <- NA
              wahl_ell6 <- out[liste6[[ende3.ell$comb.circ]], ]
              wahl_ell6[, c(11:14)] <- NA
            }else{ende3.ell <- ende2.ell}
            
            if(nrow(ende3.ell)==0){
              
              vec.6 <- which(sapply(liste6, FUN=function(X) "5" %in% X)) #
              
              ende2.ell <- ende.ell[ende.ell$comb.circ %in% vec.6, ]
              
              if(nrow(ende2.ell)!=0){
                ende3.ell <- ende2.ell[ende2.ell$sdbhd.ell==min(ende2.ell$sdbhd.ell), ][1, ]
                
                liste6 <- combn(out.vec, 6, simplify = FALSE)
                #schauen was gebraucht wird - rest NA machen
                wahl_ell <- out[5, ]
                wahl_ell[, c(11:14)] <- NA
                wahl_ell6 <- out[liste6[[ende3.ell$comb.circ]], ]
                wahl_ell6[, c(11:14)] <- NA
              }else{ende3.ell <- ende2.ell}
              
              if(nrow(ende3.ell)==0){
                
                vec.6 <- which(sapply(liste6, FUN=function(X) "6" %in% X)) #
                
                ende2.ell <- ende.ell[ende.ell$comb.circ %in% vec.6, ]
                
                if(nrow(ende2.ell)!=0){
                  ende3.ell <- ende2.ell[ende2.ell$sdbhd.ell==min(ende2.ell$sdbhd.ell), ][1, ]
                  
                  liste6 <- combn(out.vec, 6, simplify = FALSE)
                  #schauen was gebraucht wird - rest NA machen
                  wahl_ell <- out[6, ]
                  wahl_ell[, c(11:14)] <- NA
                  wahl_ell6 <- out[liste6[[ende3.ell$comb.circ]], ]
                  wahl_ell6[, c(11:14)] <- NA
                }else{ende3.ell <- ende2.ell}
                
                if(nrow(ende3.ell)==0){
                  
                  vec.6 <- which(sapply(liste6, FUN=function(X) "7" %in% X)) #
                  
                  ende2.ell <- ende.ell[ende.ell$comb.circ %in% vec.6, ]
                  
                  if(nrow(ende2.ell)!=0){
                    ende3.ell <- ende2.ell[ende2.ell$sdbhd.ell==min(ende2.ell$sdbhd.ell), ][1, ]
                    
                    liste6 <- combn(out.vec, 6, simplify = FALSE)
                    #schauen was gebraucht wird - rest NA machen
                    wahl_ell <- out[7, ]
                    wahl_ell[, c(11:14)] <- NA
                    wahl_ell6 <- out[liste6[[ende3.ell$comb.circ]], ]
                    wahl_ell6[, c(11:14)] <- NA
                  }else{ende3.ell <- ende2.ell}
                  
                  if(nrow(ende3.ell)==0){
                    
                    vec.6 <- which(sapply(liste6, FUN=function(X) "8" %in% X)) #
                    
                    ende2.ell <- ende.ell[ende.ell$comb.circ %in% vec.6, ]
                    
                    if(nrow(ende2.ell)!=0){
                      ende3.ell <- ende2.ell[ende2.ell$sdbhd.ell==min(ende2.ell$sdbhd.ell), ][1, ]
                      
                      liste6 <- combn(out.vec, 6, simplify = FALSE)
                      #schauen was gebraucht wird - rest NA machen
                      wahl_ell <- out[8, ]
                      wahl_ell[, c(11:14)] <- NA
                      wahl_ell6 <- out[liste6[[ende3.ell$comb.circ]], ]
                      wahl_ell6[, c(11:14)] <- NA
                    }else{ende3.ell <- ende2.ell}
                    
                    if(nrow(ende3.ell)==0){
                      
                      vec.6 <- which(sapply(liste6, FUN=function(X) "9" %in% X)) #
                      
                      ende2.ell <- ende.ell[ende.ell$comb.circ %in% vec.6, ]
                      
                      if(nrow(ende2.ell)!=0){
                        ende3.ell <- ende2.ell[ende2.ell$sdbhd.ell==min(ende2.ell$sdbhd.ell), ][1, ]
                        
                        liste6 <- combn(out.vec, 6, simplify = FALSE)
                        #schauen was gebraucht wird - rest NA machen
                        wahl_ell <- out[9, ]
                        wahl_ell[, c(11:14)] <- NA
                        wahl_ell6 <- out[liste6[[ende3.ell$comb.circ]], ]
                        wahl_ell6[, c(11:14)] <- NA
                      }else{ende3.ell <- ende2.ell}
                    }
                    
                  }
                  
                }
                
              }
              
            }
          }
          
        }
        
      }
    }
  }else{
    wahl_ell <- out[1, ]
    wahl_ell[1, 4:ncol(wahl_ell)] <- NA
    wahl_ell6 <- out[1:6, ]
    wahl_ell6[, 4:ncol(wahl_ell6)] <- NA
  }
  
  # wahl_circ
  # wahl_circ6
  # wahl_ell
  # wahl_ell6
  
  
  if(is.na(wahl_circ$d.circ)&is.na(wahl_ell$d.ell)){
    wahl_all_ende <- wahl_circ6[1, c("cluster", "cluster2", "cluster3")]
    wahl_all_ende$x.choosen <- NA
    wahl_all_ende$y.choosen <- NA
    wahl_all_ende$d.circ <- NA
    wahl_all_ende$d.circ2 <- NA
    wahl_all_ende$d.ell <- NA
    wahl_all_ende$d.gam <- NA
    wahl_all_ende$tegam.d_gam <- NA
    # wahl_all_ende[, colnames(wahl_circ6[19:32])] <- NA #weil die Spalten 19:32 die qgams sind
    
  }else{
    if(is.na(wahl_circ$d.circ)==F&is.na(wahl_ell$d.ell)==F){
      wahl_all_ende <- wahl_circ6[1, c("cluster", "cluster2", "cluster3")]
      hoehe <- wahl_circ6$z.unten+0.15/2 # weil Schichtbreite 0.15 ist
      model <- lm(wahl_circ6$x.circ~hoehe)
      wahl_all_ende$x.choosen <- predict(model, data.frame("hoehe"=1.3))
      model <- lm(wahl_circ6$y.circ~hoehe)
      wahl_all_ende$y.choosen <- predict(model, data.frame("hoehe"=1.3))
      model <- lm(wahl_circ6$d.circ~hoehe)
      wahl_all_ende$d.circ <- predict(model, data.frame("hoehe"=1.3))
      model <- lm(wahl_circ6$d.circ2~hoehe)
      wahl_all_ende$d.circ2 <- predict(model, data.frame("hoehe"=1.3))
      model <- lm(wahl_ell6$d.ell~hoehe)
      wahl_all_ende$d.ell <- predict(model, data.frame("hoehe"=1.3))
      model <- lm(wahl_circ6$d.gam~hoehe)
      wahl_all_ende$d.gam <- predict(model, data.frame("hoehe"=1.3))
      wahl_all_ende$tegam.d_gam <- unique(wahl_circ6$tegam.d_gam)
      # for(column in 19:32){
      #  model <- lm(wahl_circ6[, column]~hoehe)
      #  wahl_all_ende[, colnames(wahl_circ6[column])] <- predict(model, data.frame("hoehe"=1.3))
      # }
    }
    
    if(is.na(wahl_circ$d.circ)==F&is.na(wahl_ell$d.ell)){
      wahl_all_ende <- wahl_circ6[1, c("cluster", "cluster2", "cluster3")]
      hoehe <- wahl_circ6$z.unten+0.15/2 # weil Schichtbreite 0.15 ist
      model <- lm(wahl_circ6$x.circ~hoehe)
      wahl_all_ende$x.choosen <- predict(model, data.frame("hoehe"=1.3))
      model <- lm(wahl_circ6$y.circ~hoehe)
      wahl_all_ende$y.choosen <- predict(model, data.frame("hoehe"=1.3))
      model <- lm(wahl_circ6$d.circ~hoehe)
      wahl_all_ende$d.circ <- predict(model, data.frame("hoehe"=1.3))
      model <- lm(wahl_circ6$d.circ2~hoehe)
      wahl_all_ende$d.circ2 <- predict(model, data.frame("hoehe"=1.3))
      wahl_all_ende$d.ell <- NA
      model <- lm(wahl_circ6$d.gam~hoehe)
      wahl_all_ende$d.gam <- predict(model, data.frame("hoehe"=1.3))
      wahl_all_ende$tegam.d_gam <- unique(wahl_circ6$tegam.d_gam)
      # for(column in 19:32){
      #  model <- lm(wahl_circ6[, column]~hoehe)
      #  wahl_all_ende[, colnames(wahl_circ6[column])] <- predict(model, data.frame("hoehe"=1.3))
      # }
      
    }
    if(is.na(wahl_circ$d.circ)&is.na(wahl_ell$d.ell)==F){
      wahl_all_ende <- wahl_circ6[1, c("cluster", "cluster2", "cluster3")]
      hoehe <- wahl_ell6$z.unten+0.15/2 # weil Schichtbreite 0.15 ist
      model <- lm(wahl_ell6$x.ell~hoehe)
      wahl_all_ende$x.choosen <- predict(model, data.frame("hoehe"=1.3))
      model <- lm(wahl_ell6$y.ell~hoehe)
      wahl_all_ende$y.choosen <- predict(model, data.frame("hoehe"=1.3))
      wahl_all_ende$d.circ <- NA
      wahl_all_ende$d.circ2 <- NA #kann man ueberlegen, ob man den circ nicht hier hineintut - denn es gibt ihn ja auch wenn es den circ nicht gibt
      model <- lm(wahl_ell6$d.ell~hoehe)
      wahl_all_ende$d.ell <- predict(model, data.frame("hoehe"=1.3))
      model <- lm(wahl_ell6$d.gam~hoehe)
      wahl_all_ende$d.gam <- predict(model, data.frame("hoehe"=1.3))
      wahl_all_ende$tegam.d_gam <- unique(wahl_ell6$tegam.d_gam)
      # for(column in 19:32){
      #  model <- lm(wahl_ell6[, column]~hoehe)
      #  wahl_all_ende[, colnames(wahl_ell6[column])] <- predict(model, data.frame("hoehe"=1.3))
      # }
    }
  }
  
  
  return(wahl_all_ende)
}

fineCluster <- function(fileFinder, dbhPath, allDBHs = FALSE, nr_cores = 0, 
                        referenced = FALSE, bushPreparation = FALSE, filterSOR = FALSE,
                        cutWindow = c(-1000, -1000, 2000), silent = FALSE, dirPath = paste0(getwd(), "/")){
  cat("\nStarting fineCluster(), ")
  fineTime1 <- Sys.time()
  cat("Current time is", format(fineTime1, "%H:%M:%S"), "\n")
  groundPath <- v.env$groundPath
  
  #cat("Handling set", file, "-", file.name.i, "\n")
  end_path <- paste0(dbhPath, "fineCluster/")
  file.list <- list.files(end_path)
  
  cat(" -> analyzing", length(file.list), "clusters...\n")
  
  
  if(nr_cores == 0){
    #cat("Automatic number of core setting:\n   ")
    nr_cores <- detectCores()
    #cat(nr_cores, "cores available - setting to max ")
    nr_cores <- nr_cores - 1
    
    # maxCores <- 10
    # cat(maxCores, "cores recommended.\n")
    # if(maxCores < 1){
    #   maxCores <- 1
    # }
    # if(maxCores <=  nr_cores){
    #   nr_cores <- maxCores
    # }
  }
  
  
  if(nr_cores == 1){
    
    cat("Going for serial cluster analysis:")
    for(datei in 1:length(file.list)){
      if(datei%%20 == 1) cat("\n    ")
      cat(paste0("", datei, "-"))
      nowCircles <- fineCluster_i(circleFile = paste0(end_path, file.list[datei]))
      
      if(datei==1){
        #entdeckung <- wahl
        entdeckung_all <- nowCircles
      }else{
        #entdeckung <- rbind(entdeckung, wahl)
        entdeckung_all <- rbind(entdeckung_all, nowCircles)
      }
    }
  }
  else 
  {
    # PARALLEL PART
    cat(" -> going now parallel on", nr_cores, "cores...")
    
    library(doParallel)
    registerDoParallel(nr_cores)
    
    file_parallelProtocol <- paste0(dbhPath, "temp_par_fineClust.txt")
    file.create(file_parallelProtocol)
    fdc <<- foreach(datei = 1:length(file.list),  .errorhandling = 'remove', 
                    .export='fineCluster_i') %dopar% {
                      sink(file_parallelProtocol, append = T)
                      cat(paste0(datei, "-"))
                      nowCircles <- fineCluster_i(circleFile = paste0(end_path, file.list[datei]))
                      
                      cat(paste0("-", datei, "\n"))
                      sink()
                      return(nowCircles)
                    }
    
    stopImplicitCluster()
    
    unlink(file_parallelProtocol)
    
    skippedTrees <- sum(lengths(fdc) == 0)
    if(skippedTrees > 0){
      cat("There were", skippedTrees, "clusters removed.\n")
      
      table(lengths(fdc))
      fdc2 <- NULL
      for(i in 1:length(fdc)){
        if(!is.null(fdc[[i]])){
          if(is.null(fdc2)){
            fdc2 <- fdc[[i]]
          } else {
            fdc2 <- rbind(fdc2, fdc[[i]])
          }
        }
      }
      df <- fdc2
      cat("There were", length(df[,1]), "trees calculated on our point.")
    } else {
      df <- data.frame(matrix(unlist(fdc), nrow=length(fdc), byrow=TRUE))
      colnames(df) <- colnames(fdc[[1]])
      for(j in c(1:length(df[1,]))){ # all num 
        df[,j] <- as.numeric(as.character(df[,j]))
      }
      # #special case, convert treeName if it works
      # tempTreeNames <- as.numeric(as.character(df[,2]))
      # if(sum(is.na(tempTreeNames))==0){
      #   df[,2] <- tempTreeNames
      # }
    }
    entdeckung_all <- df
  }
  
  
  #gleich filtern auf die Baeme die es auch tats. gibt
  entdeckung_all <- entdeckung_all[is.na(entdeckung_all$x.choosen)==F, ]
  write.csv2(entdeckung_all, paste(dbhPath, "clusters_raw.csv", sep=""), row.names = F)
  
  #write.csv2(entdeckung_all, paste("F:/Testen_Intensitaet/END_TLS/", "254", ".csv", sep=""), row.names = F)
  
  
  
  
  #} # only one file
  
  
  
  cat("done!\n")
  cat("\nCreating output tree cluster list:\n")
  ### want to create output file with especially z coordinates, then want to leave referencing for later (when all is done)
  
  #  setStr <- generateSetString(fileFinder = fileFinder, mode = mode,
  #                clipHeight = clipHeight, bottomCut = bottomCut,
  #                bushPreparation = bushPreparation,
  #                filterSOR = filterSOR, cutWindow = cutWindow, silent = TRUE)
  #  dbhPath <- paste0(dirPath, setStr, "_dbh/")
  # trees.file <- paste0(dbhPath, "clusters_raw.csv")
  # cat("Reading detected stems from ", trees.file, "\n", sep = "")
  # trees <- read.csv2(trees.file)
  trees <- entdeckung_all
  trees <- trees[order(trees$d.gam, decreasing = TRUE), ]
  
  
  ## STRIPING ESSENTIAL COLUMNS to get x and y column
  {
    if(sum(colnames(trees)=="x")==0){
      # there is no column of x, let's find an alternative:
      colnames(trees)[colnames(trees)=="x.choosen"] <- "x"
    }
    if(sum(colnames(trees)=="y")==0){
      # there is no column of y
      colnames(trees)[colnames(trees)=="y.choosen"] <- "y"
    }
    if(sum(colnames(trees)=="id")==0){
      # there is no column of id
      trees$id <- c(1:length(trees[, 1]))
      cat(" * assigning new unique id from 1 to", length(trees[, 1]), "\n")
      #colnames(trees)[colnames(trees)=="cluster"] <- "id"
    }
  }
  
  if(sum(colnames(trees)=="x") + sum(colnames(trees)=="y") + sum(colnames(trees)=="id") != 3){
    cat("Detected tree file is missing an input value (either x, y or id!\n")
    cat("Terminating alloation, please check your files.\n\n")
    return()
  }
  
  cat(" * attaching z-values from reading _ground_min.grd model\n")
  tryCatch(
    {
      # read in raster file
      dtm_z <- raster(paste0(dirPath, groundPath, fileFinder, "_ground_min.grd"))
    }, error = function(error_condition) {
      cat("Error in reading the min dtm-model!")
      return()
    })
  trees$z <- round(extract(dtm_z, SpatialPoints(data.frame(x = trees$x, y = trees$y))) + 1.3, 3)
  #extracting the z-value from a dtm, was hard work to find that out... so easy
  #    extract(dtm_c, SpatialPoints(data.frame(x=15, y=-10)))
  
  
  
  cat(" * drawing new random stem ids\n")
  rnlist <- sample(1:65535, length(trees[, 1]))
  
  
  mergedSet <- data.frame("x" = round(trees$x, 3), "y" = round(trees$y, 3),
                          "z" = trees$z, "id" = trees$id)
  if(sum(colnames(trees) == "cluster")==1){
    mergedSet$cluster <- trees$cluster
  }
  mergedSet$randomCol <- rnlist
  mergedSet$dbh <- round(trees$d.gam, 1)
  if(allDBHs){
    mergedSet$d.circ <- round(trees$d.circ, 3)
    mergedSet$d.circ2 <- round(trees$d.circ2, 3)
    mergedSet$d.ell <- round(trees$d.ell, 3)
    mergedSet$d.gam <- round(trees$d.gam, 3)
    mergedSet$tegam.d_gam <- round(trees$tegam.d_gam, 3)
  }
  
  write.table(mergedSet, file = paste0(dbhPath, "trees_dbh.txt"),
              row.names = FALSE, sep = "\t")
  
  cat(paste0("\nFile ", "\"", dbhPath, "trees_dbh.txt", "\""), "created.\n")
  
  
  fineTime2 <- Sys.time()
  cat("Fine clustering done in a ")
  print.difftime(round(fineTime2 - fineTime1, 1))
  cat("\n")
  
}



