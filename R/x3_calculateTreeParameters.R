
measureNewDBH <- function(treeLAS, nowMeta, picPath = getwd(),
                          threshold_percent = 50,
                          threshold_percent_kde = 94){

  library("lidR")
  library("ks")
  library("mgcv")
  library("stats")
  library("treeX")

  if(1==2){
    threshold_percent = 50
    threshold_percent_kde = 94
  }



  dbhMeta <- data.frame("id" = nowMeta$id, "dbh_ref" = nowMeta$dbh)
  dbhMeta$dbh <- -1
  dbhMeta$dbh.gam <- -1
  dbhMeta$sd.gam <- -1
  dbhMeta$var.dbh <- -1
  dbhMeta$var.dbh.gam <- -1
  dbhMeta$var.sd.gam <- -1
  dbhMeta$x <- -1
  dbhMeta$x.ref <- nowMeta$x
  dbhMeta$y <- -1
  dbhMeta$y.ref <- nowMeta$y
  dbhMeta$slices <- -1
  dbhMeta$move <- -1


  cat("Remeasure DBH properly - intensity filtering",threshold_percent,"% and KDE filtering per slice",threshold_percent_kde,"%.\n")

  # clipping to z-value
  treeLAS <- filter_poi(treeLAS, Z > nowMeta$z - 0.3, Z < nowMeta$z + 0.3)
  # suffice: 6 slices from 1 m to 1.60 m
  nowId <- nowMeta$id

  # INTENSITY FILTERING
  threshold <- quantile(treeLAS@data$Intensity, 1 - threshold_percent / 100) # all above are 5 %
  clustInt <- filter_poi(treeLAS, Intensity > threshold)
  #cat(" ", clustInt@header@PHB$`Number of point records`, "pts ")

  # WEIGHED KERNEL DENSITY 2D
  #image(kde$eval.points[[1]], kde$eval.points[[2]], kde$estimate,
  # col = viridis::viridis(20), xlab = "x", ylab = "y")
  kde_dat <- data.frame("x" = clustInt@data$X, "y" = clustInt@data$Y)
  if(length(kde_dat$x) < 10){
    cat("\n",dbhMeta$id, "- too few points in seed, skip this one\n")
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




  thisId <- ifelse(is.null(nowMeta$treeName), nowMeta$id, nowMeta$treeName)

  par(mfrow = c(1, 1))
  {
    png(filename = paste0(picPath, thisId, "_circles.png"), width = 1400, height = 800)
    par(mfrow = c(3, 4), oma = c(1,1,4,1))
    circles <- NULL

    cat("slice ")
    dbhStart <- Sys.time()
    for(j in 1:length(slices)){
      cat(j, "- ")
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
        points(dbhMeta$x.ref, dbhMeta$y.ref, col = "grey", cex = 2, pch = 4)

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

          mtext(paste0(" move ", round(100*sqrt((x.circle - dbhMeta$x.ref)^2 + (y.circle - dbhMeta$y.ref)^2), 1)),
                side = 1, adj = 0, line = -2, col = "grey", cex = 1.5)
          draw.circle(x.circle, y.circle, dbh.circle / 200, border = "red")

          points(x.circle, y.circle, col = "red", cex = 2, pch = 4)
          lines(x = c(x.circle, dbhMeta$x.ref),
                y = c(y.circle, dbhMeta$y.ref), col = "grey")

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
      paste0(dbhMeta$id,
             " dbhs: circ ", round(median(circles[3, ])*200, 1),
             " (vr", round(sd(circles[3, ]*200, na.rm = T),2), ")",
             " gam ", round(median(circles[6, ], na.rm = T), 1),
             " (vr", round(sd(circles[6, ], na.rm = T),2), ")",
             " sd ", round(median(circles[7, ], na.rm = T), 1),
             " (vr", round(sd(circles[7, ], na.rm = T),2), ")",
             " --- alt ", round(dbhMeta$dbh_ref, 1))
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

      dbhMeta$dbh <- dbh.mean
      dbhMeta$dbh.gam <- dbh.gam.mean
      dbhMeta$sd.gam <- gam.sd.mean
      dbhMeta$var.dbh <- round(sd(circles[3, ]*200, na.rm = T),2)
      dbhMeta$var.dbh.gam <- round(sd(circles[6, ], na.rm = T),2)
      dbhMeta$var.sd.gam <- round(sd(circles[7, ], na.rm = T),4)
      # HUOM: Not real variance, but sd just to prevent name confusion with gam sd (scattering)
      dbhMeta$x <- x.mean
      dbhMeta$y <- y.mean
      dbhMeta$slices <- length(slices)
      dbhMeta$move <- round(100*sqrt((x.mean - dbhMeta$x.ref)^2 + (y.mean - dbhMeta$y.ref)^2), 1)
    } else {
      #timeNB <- as.difftime(dbhStop - dbhStart)
      #cat("No DBH estimated in", round(timeNB, 1), units(timeNB), "\n")
    }

    circles <- rbind(circles, slices)
    circles <- rbind(circles, slices + 0.1)
    rownames(circles) <- c("cl.x", "cl.y", "cl.d", "cl.ka1", "cl.ka2", "cl.dgam", "cl.sdgam", "lower", "upper")

    #dev.off()

  }
  #selectorius <- c(1:4, 8, 10, 12, 14)
  #cat(paste(dbhMeta[selectorius]), sep = "\t")
  #cat("\n")
  write.table(dbhMeta, paste0(picPath, nowId, "_rem_dbhs.txt"), row.names = F, sep = "\t", dec = ".")
  write.csv2(t(circles), paste0(picPath, nowId, "_circles_x6.csv"), row.names = F)



  {
    #cat("\nDeciding now, which new diameters to take as valid.\n")
    # check which new measured to take into
    # LADIWALDI RULEBOOK:
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

    allTrees <- T
    if(allTrees){
      improveList <- dbhMeta
    } else {
      improveList <- clustList[clustList$id >= 9000,]
    }
    improveList$dbh.gam[improveList$dbh.gam > 140] <- -11
    improveList$dbh.gam[improveList$dbh.gam < 0] <- -11
    improveList$dbh.gam[is.na(improveList$dbh.gam)] <- -11
    improveList$dbh[improveList$dbh < 0] <- -11
    improveList$dbh[is.na(improveList$dbh)] <- -11
    #cat("In total there are", length(improveList$id),
    #    "new measured dbhs ( na.gam =",sum(improveList$dbh.gam == -11),
    #    "and na.circ =", sum(improveList$dbh == -11),").\n")

    improveList$diff.dist <- sqrt((improveList$x.ref - improveList$x)^2 + (improveList$y.ref - improveList$y)^2 )
    improveList$diff.dbh <- abs(improveList$dbh_ref - improveList$dbh)
    improveList$diff.dbh.gam <- abs(improveList$dbh_ref - improveList$dbh.gam)

    #cat("5 cm accepted tolerance to guessed DBH: ")
    acceptToleranceCM <- 5
    ##table(improveList$dbh.gam != -11)
    improveList$dbh.check <- "Fi" #good ones
    improveList$dbh.sugg <- improveList$dbh.gam
    improveList$dbh.bad <- " "
    improveList$dbh.bad[improveList$diff.dbh.gam > acceptToleranceCM] <- "X"
    ##improveList$dbh.bad[improveList$dbh.gam == -11] <- "X"
    #table(improveList$dbh.bad[improveList$dbh.gam != -11])

    #cat(sum(improveList$dbh.bad[improveList$dbh.gam != -11] != "X"), "assigned\n")
    #table(improveList$dbh.bad)



    #cat("90-110 % accepted ratio to guessed DBH: ")
    ##set 90 - 110 % circle/gam equivalence of both measures as acceptable (Buche)
    improveList$dbh.check[improveList$dbh.bad == "X"] <- "Bu" #65 with 90-110% equivalence
    improveList$dbh.est.ratio <- improveList$dbh/improveList$dbh.gam*100
    improveList$dbh.bad[improveList$dbh.est.ratio >= 90 &
                          improveList$dbh.est.ratio <= 110 &
                          improveList$dbh.gam != -11] <- " "
    #table(improveList$dbh.bad[improveList$dbh.gam != -11])
    #cat(sum(improveList$dbh.bad[improveList$dbh.gam != -11] != "X"), "assigned\n")
    #table(improveList$dbh.bad)




    #cat("5 cm allowed discrepancy to circle DBH: ")
    ##set absolute difference 5 cm as acceptable (Ulme)
    improveList$dbh.check[improveList$dbh.bad == "X"] <- "Ul" #6 with 5cm difference
    improveList$dbh.est.ratio.diff <- abs(improveList$dbh - improveList$dbh.gam)
    improveList$dbh.bad[improveList$dbh.est.ratio.diff <= acceptToleranceCM &
                          improveList$dbh.gam != -11] <- " "
    #table(improveList$dbh.bad[improveList$dbh.gam != -11])
    ##a_Rosalia: 1 bad ones + 133
    ##b_Altmuenster: 2 bad ones + 58
    #cat(sum(improveList$dbh.bad[improveList$dbh.gam != -11] != "X"), "assigned\n")
    #table(improveList$dbh.bad)



    #cat("CIRCLE 5 cm discrepancy to guessed DBH: ")
    ## 133 bad gam ones: take the dbh of circular measurements if difference is less than 5 cm
    improveList$dbh.check[improveList$dbh.bad == "X"] <- "Ta" #51 dbh circle values
    choice.bad.takeCircledbh <- improveList$dbh.bad == "X" & improveList$diff.dbh <= 5
    improveList$dbh.bad[choice.bad.takeCircledbh] <- " "
    #table(improveList$dbh.bad) #82 left
    improveList$dbh.sugg[choice.bad.takeCircledbh] <- improveList$dbh[choice.bad.takeCircledbh]
    #table(improveList$dbh.bad) #82 no good measurement, Altmuenster: 18 no good
    #cat(sum(improveList$dbh.bad != "X"), "assigned\n")
    #improveList[improveList$dbh.bad=="X",]




    # trees smaller than 1 cm dbh
    #
    thoseSmallerThanOneCm <- improveList$dbh.sugg != -11 & improveList$dbh.sugg < 1
    if(sum(thoseSmallerThanOneCm) > 0){
      #cat("Trees smaller than 1 cm correction for",sum(thoseSmallerThanOneCm),"trees.")
      improveList[improveList$dbh.sugg != -11 & improveList$dbh.sugg < 1,]
      improveList$dbh.sugg[improveList$dbh.sugg != -11 & improveList$dbh.sugg < 1] <- 1 # set to 1 cm

    }
    improveList$dbh.check[improveList$dbh.sugg == -11] <- "Ks" #82 not measured, stay estimates
    #print(table(improveList$dbh.check))
    #sum(improveList$dbh.gam == -11) #133 no gam
    #sum(improveList$dbh.sugg == -11) #82 not measured (51 from circle dbh)
    #notImprovedList <- improveList[improveList$dbh.sugg == -11,]
    ##table(notImprovedList$dbh_ref) # only up to 5 cm dbh, Altmuenster one very outside 30 cm (not anymore in area)



    #cat("There are",length(notImprovedList$dbh),"new trees without improved dbh.\n")
    #cat(notImprovedList$id, sep = " - ")
    #cat("\n")


    ##improveList <- improveList[improveList$dbh.gam != -1,]
    ##improveList$id[improveList$dbh.check == "X"] # check these ones, they have
    ##improveList[improveList$dbh.sugg > 100,]





    improveList$shift.x <- NA
    improveList$shift.y <- NA

    # SHIFT ALL
    omitShifting <- TRUE
    if(!omitShifting){
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

      #
      # dbhMeta$x[dbhMeta$x == -1] <- NA
      # dbhMeta$y[dbhMeta$y == -1] <- NA

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
    } else {

      mergeImprove <- data.frame("id" = improveList$id,
                                 "dbh.ref" = improveList$dbh_ref,
                                 "dbh.sugg" = improveList$dbh.sugg,
                                 "dbh.check" = improveList$dbh.check)
    }

  }



  if(mergeImprove$dbh.sugg == -11){
    returnDF <- data.frame("x" = 0,
                           "y" = 0,
                           "dbh" = mergeImprove$dbh.ref)
    return(returnDF)
  } else {
    cat("Changing dbh from", mergeImprove$dbh.ref, "to",
        round(mergeImprove$dbh.sugg, 2),"cm!\n")
    returnDF <- data.frame("x" = dbhMeta$x,
                           "y" = dbhMeta$y,
                           "dbh" = mergeImprove$dbh.sugg)
    return(returnDF)
  }

}





#' calculates the individual tree parameters for dedicated tree
#'
#' @param tree number of StemID which to measure
#' @param treeLAS.path individual Tree LAS file
#' @param treeName additional name for the images (specific trees)
#' @param z_stemBase predefined z-coordinate for start of the stem, overrides nowMeta$z
#' @param nowMeta additional information about tree species, xyz coordinates and dbh
#' @param detail.level how many points shall be kept of crown file for hulling 0 = standard (2 pts every 25x25cm), -1 = rough (5 pts per 1x1m), -3 very rough (5 pts per 3x3m), 3 very fine (3 pts every 10x10cm)
#' TODO: implement the crownBase list input with passing the column crownBase in nowMeta, I think I fixed it in 23-11 (without testing)
#' @export
computeTree_i <- function(treeLAS.path,
                          treeName = "",
                          drawImage = F,
                          z_stemBase = NA,
                          nowMeta = NA,
                          detail.level = 0,
                          decimateTreeForFasterCrowns = FALSE,

                          #writeLAZ = TRUE, writePicture = TRUE,
                          #density.CrownLAS = 800, # how many points in random() filter for LAS to be kept for crown hulling

                          vol.alpha = 2, area.alpha = 0.3,
                          alternativeCrownBase.Ratio = 0.3, fogFilter.estHeight = FALSE,
                          measurePath = dirname(treeLAS.path),
                          fileFinder = "",

                          remeasureDBH = FALSE,
                          do.Plot = TRUE,
                          limitSpanSide = 30, limitSpanArea = 440,
                          referenceDiameterLimit = 150
                          ){

  # cat("Single Tree Measurements in folder\n   ", measurePath,"\n")
  gstart <- Sys.time()

  {

    exceededReferenceLimit <- 0

    crownImagePath <- paste0(measurePath,"/paramPics/")
    stemAnalysisPath <- paste0(measurePath,"/stemAnalysis/")
    if(!dir.exists(crownImagePath)) dir.create(crownImagePath)
    if(!dir.exists(stemAnalysisPath)) dir.create(stemAnalysisPath)


    cat("Reading file", basename(treeLAS.path), "...")
    if(!file.exists(treeLAS.path)){
      cat("File not found! Skipping...\n\n")
      return(data.frame())
    }
    treeLAS <- readLAS(treeLAS.path)
    cat(thMk(treeLAS@header$`Number of point records`), "pts ok.\n")





    if(treeName == ""){
      #treeName <- paste0(sprintf("%04d", trees[i]),".laz")
      treeName <- strRep(basename(treeLAS.path), ".las", "")
      treeName <- strRep(treeName, ".laz", "")
    }


    metaVars <- data.frame("file"=treeLAS.path,
                           "treeName"=treeName,
                           "nPoints" = treeLAS@header@PHB$`Number of point records`,
                           "est.z.StemBase" = nowMeta$z-1.3,
                           "est.DBH"=0, "est.height"=0,
                           "est.crownLength"=0, "est.crownBase"=0,
                           "est.crownArea"=0, "est.crownAreaAlt"=0,
                           "est.crownVolume"=0,
                           #"est.crownDiameter"=0, "est.DiamCrownMax"=0, "est.DiamCrownMax.Angle"=0,
                           #"est.DiamCrownMin"=0, "est.DiamCrownMin.Angle"=0,
                           "est.x.DBH"=0, "est.y.DBH"=0, "est.z.DBH"=0,
                           #"est.x.crownDiameter"=0, "est.y.crownDiameter"=0,
                           stringsAsFactors=FALSE) #not implemented yet, was interesting!

    if(remeasureDBH){
      picPath <- paste0(measurePath, "/stemCircles/")
      if(!dir.exists(picPath)){
        dir.create(picPath)
      }
      nowMeta$treeName <- treeName
      measureDF <- measureNewDBH(treeLAS, nowMeta, picPath = picPath)
      if(measureDF$dbh < 100 && measureDF$dbh > 0 ){
        metaVars$est.DBH <- measureDF$dbh
        nowMeta$dbh <- measureDF$dbh
        if(!(measureDF$x == 0 && measureDF$y == 0)){
          metaVars$est.x.DBH <- measureDF$x
          metaVars$est.y.DBH <- measureDF$y
          nowMeta$x <- measureDF$x
          nowMeta$y <- measureDF$y
        }
      }
    }

    if(drawImage){
      up.h <- treeLAS@header@PHB$`Max Z`
      if(fogFilter.estHeight){
        up.h <- quantile(treeLAS@data$Z, 0.9999, names = F)
      }
      down.h <- nowMeta$z - 1.3
      diffDraw <- treeLAS@header@PHB$`Max Z` - up.h
      maxPoint <- treeLAS@data[which(treeLAS@data$Z == treeLAS@header@PHB$`Max Z`),c(1:3)]
      maxPoint <- maxPoint[1,]
      metaVars$est.height <- round(up.h - down.h, 2)

      metaVars$topCrownX <- maxPoint$X
      metaVars$topCrownY <- maxPoint$Y
      metaVars$topCrownZ <- maxPoint$Z

      metaVars$topLeanDist <- sqrt((maxPoint$Y - nowMeta$y)^2 + (maxPoint$X - nowMeta$x)^2)
      metaVars$topLeanAngle <- atan2((maxPoint$Y - nowMeta$y), (maxPoint$X - nowMeta$x))/2/pi*360

      crownPathSingles <- paste0(dirname(treeLAS.path), "/")
      if(metaVars$nPoints > 1000 & decimateTreeForFasterCrowns){
        # if(smallTrees){
        #   treeLAS <- decimate_points(treeLAS, random(10000))
        # } else {
          treeLAS <- decimate_points(treeLAS, random(4000))

      }
      # Picture of tree from three sides with height measurement
      png(paste0(crownPathSingles, treeName,".png"), height = 1024, width = 2048, type = "cairo")
      tryCatch(
        {
          par(mfrow=c(1, 3), oma = c(0, 0, 3, 0))
          plot(treeLAS@data$X, treeLAS@data$Z, asp = 1, cex = 0.001, ylim = c(down.h, up.h))
          abline(h = c(up.h, down.h), col = "red")
          lines(x = c(nowMeta$x, maxPoint$X), y = c(down.h, up.h), col = "green", lty = 2)
          plot(treeLAS@data$Y, treeLAS@data$Z, asp = 1, cex = 0.001, ylim = c(down.h, up.h))
          abline(h = c(up.h, down.h), col = "red")
          lines(x = c(nowMeta$y, maxPoint$Y), y = c(down.h, up.h), col = "green", lty = 2)
          plot(treeLAS@data$X, treeLAS@data$Y, asp = 1, cex = 0.001)
          lines(x = c(nowMeta$x, maxPoint$X), y = c(nowMeta$y, maxPoint$Y), col = "green", lwd = 2)
          plotTitle.sin <- paste0(fileFinder, "   ", treeName, "   dbh = ", round(nowMeta$dbh,2), "cm   est.h = ",metaVars$est.height, "m")
          if(is.element("species", colnames(nowMeta))){
            if(!is.null(nowMeta$species)){
              plotTitle.sin <- paste0(fileFinder, "   ", nowMeta$species,"   ", treeName, "   dbh = ", round(nowMeta$dbh,2), "cm   est.h = ",metaVars$est.height, "m")
            }
          } else if(is.element("Bart", colnames(nowMeta))){
            if(!is.null(nowMeta$Bart)){
              plotTitle.sin <- paste0(fileFinder, "   ", nowMeta$Bart,"   ", treeName, "   dbh = ", round(nowMeta$dbh,2), "cm   est.h = ",metaVars$est.height, "m")
            }
          }

          mtext(plotTitle.sin, outer = TRUE, cex = 2.0)
          drawNumberPoints <- T
          if(drawNumberPoints){
            mtext(paste0("nPoints = ", thMk(metaVars$nPoints)), outer = TRUE, cex = 1.6, line = -2.7)
          }

        }, error=function(e){
          cat("Error when writing the image, closing device.\n")
        })
      dev.off()
      #par(mfrow=c(1,1))

    }



    if(is.null(nrow(nowMeta))){
      cat("No meta-data specified.\n")
      nowMeta <- data.frame("file" = treeLAS.path, "treeName" = treeName)
    }

    if(!is.na(z_stemBase)){
      cat("Stem base at z =", round(z_stemBase, 1), "m (z_stemBase)\n")
      nowMeta$z <- z_stemBase
    } else {
      if(is.element("z", colnames(nowMeta))){
        cat("Stem base at z =", round(nowMeta$z[1], 1)-1.3, "m (nowMeta$z)\n")
      } else if(is.element("Z", colnames(nowMeta))) {
        cat("Stem base at z =", round(nowMeta$Z[1], 1)-1.3, "m (nowMeta$Z)\n")
        nowMeta$z <- nowMeta$Z
      } else {
        nowMeta$z <- 0
        nowMeta$z <- quantile(treeLAS@data$Z, 0.01) + 1.3
        cat("Stem base at z =", round(nowMeta$z[1], 1) - 1.3, "m (estimated)\n")
        cat("  -> we used 1 % quantile of z-coordinates.\n")
      }
    }

  }

  if(is.element("crownBase", colnames(nowMeta))){
    cat("Manual Crown base at z =", round(nowMeta$crownBase, 2), ".\n")
  }


  # makes crown hulling more efficient, for really big tree we reduce it even more
  #treeTooBig <- FALSE

  tryCatch({




    # calculate crownParameters
    {

      isTree <- TRUE #might come useful to identify not tree files, not implemented yet
      #plot(treeLAS)

      # removing the gps coordinates very tall numbers
      coordXtop <- treeLAS@header@PHB$`Max X`
      if(coordXtop > 1000){
        shX <- trunc(coordXtop/1000)*1000
        treeLAS@data$X <- treeLAS@data$X - shX
      }
      coordYtop <- treeLAS@header@PHB$`Max Y`
      if(coordYtop > 1000){
        shY <- trunc(coordYtop/1000)*1000
        treeLAS@data$Y <- treeLAS@data$Y - shY
      }
      if(coordYtop > 1000 || coordXtop > 1000){
        treeLAS <- LAS(data = treeLAS@data)
      }
      #treeLAS@data$Z <- treeLAS@data$Z # +100 to get all positive points.




      #### SETTINGS FOR CROWN BASE DETECTION ####
      # the criteria for crown base is 1m of more then 2x reference diameter
      filter.quant <- 0.15 # removing 15 % of all points too far away from center x and y
      tolerance <- 0.1 #adding 10 cm again after quantile filtering not to cut off everything...
      slice.height <- 0.05 #in a 5 cm slice every diameter is searched to find crown base upwards

      if(fogFilter.estHeight){
        cat("Fog and snow filter applied to height estimation - working with 99,99 % quantile instead of max Z!\n")
      }


      if(is.element("crownBase", colnames(nowMeta))){
        cat("No calculation of crown bases - we use input from nowMeta:  ")
        cat(nowMeta$crownBase, "\n")
        # Reading in crown base file
        #crownBaseList <- read.table(path.inputCrownBase, sep = "\t", header = TRUE)
        #cat(" - containing",nrow(crownBaseList),"file entries.\n\n")

      } else {
        #cat("Crown bases will be calculated!\n")
        #alternativeCrownBase.Ratio <- 0.3 # of tree height, set above
        #cat("If no crown base is found, the alternative measures assume", alternativeCrownBase.Ratio*100, "% of tree height.\n")

      }






      # initializes
      # minZ <- min(treeLAS@data$Z)
      minZ <- nowMeta$z-1.3



      # ESTIMATE HEIGHT

      cat("height is ")
      up.h <- treeLAS@header@PHB$`Max Z`
      if(fogFilter.estHeight){
        cat(" - snow reduced -")
        up.h <- quantile(treeLAS@data$Z, 0.9999, names = F)
      }

      down.h <- nowMeta$z-1.3
      maxPoint <- treeLAS@data[which(treeLAS@data$Z == treeLAS@header@PHB$`Max Z`),c(1:3)]
      maxPoint <- maxPoint[1,]
      metaVars$est.height <- round(up.h - down.h, 2)

      altCrownBase <- (treeLAS@header@PHB$`Max Z` - minZ) * alternativeCrownBase.Ratio
      cat(metaVars$est.height, "m (z.alt.cr.b =", round(altCrownBase, 3), "m)\n")
      altCrownBase <- minZ + altCrownBase





      # we use these heights to get the mean stem axis point
      heightGrip <- seq(from = 0.0, to = max(treeLAS@data$Z) - minZ, by = slice.height)
      centers <- data.frame("x" = nowMeta$x, "y" = nowMeta$y, "z" = minZ + 1.3,
                            "d" = nowMeta$dbh/100, "lower" = 0)
      firstCenter <- centers
      stems <- data.frame("x" = nowMeta$x, "y" = nowMeta$y, "z" = minZ,
                          "d" = nowMeta$dbh/100, "movedCM" = 0, "angle" = 0, "lower" = 0)

      crownCounter <- 0
      groundCutHeight <- 0 # if floor was removed, set dbh higher




      tempTree <- decimate_points(treeLAS, random(1000))
      #plot(tempTree)



      {
        plot.now <- FALSE
        t1 <- Sys.time()
        altCrownLAS <- filter_poi(treeLAS, Z >= altCrownBase)
        if(plot.now) tempTree <- altCrownLAS



        # # possibility to do noise filter
        # altCrownLAS <- classify_noise(altCrownLAS, ivf(res = 0.2, n = 100))
        # plot(altCrownLAS, color = "Classification")
        # altCrownLAS <- filter_poi(altCrownLAS, Classification != LASNOISE)

        altCrownLAS@data$Z <- 0
        #altCrownLAS <- decimate_points(altCrownLAS, random(density.CrownLAS))
        #plot(altCrownLAS@data$Y ~altCrownLAS@data$X, cex = 0.0001, asp = 1)
        #altCrownLAS <- voxelize_points(altCrownLAS, 0.3)
        #altCrownLAS <- decimate_points(altCrownLAS, random_per_voxel(res = 0.10))
        altCrownLAS <- decimate_points(altCrownLAS, homogenize(density = 30, res = 0.2))
        #altCrownLAS <- decimate_points(altCrownLAS, random(density.CrownLAS))
        altCrownLAS <- filter_duplicates(altCrownLAS)



        #crownOutlineXYALT <- ahull(xysetALT$x, xysetALT$y, alpha = area.alpha)
        crownOutlineXYALT <- ahull(altCrownLAS@data$X, altCrownLAS@data$Y, alpha = area.alpha)
        t2 <- Sys.time()
        #print.difftime(t2-t1)

        if(plot.now)  plot(tempTree@data$Y ~ tempTree@data$X, cex = 0.0001, asp = 1)
        if(plot.now)  points(altCrownLAS@data$Y ~ altCrownLAS@data$X, cex = 0.5, pch = "+", asp = 1, col = "red")
        if(plot.now)  plot(crownOutlineXYALT, add = T, col = "red", wpoints = F)

      }


      # crownOutlineXYALT <- ahull(xysetALT$x, xysetALT$y, alpha = 0.3)
      metaVars$est.crownAreaAlt <- round(areaahulleval(crownOutlineXYALT),1)




      spanX <- treeLAS@header@PHB$`Max X` - treeLAS@header@PHB$`Min X`
      spanY <- treeLAS@header@PHB$`Max Y` - treeLAS@header@PHB$`Min Y`
      if(spanX > limitSpanSide || spanY > limitSpanSide || spanX * spanY > limitSpanArea){
        #exceededReferenceLimit <- exceededReferenceLimit + 1
        #treeTooBig <- T
        cat("Very big lateral span for", treeName, " - no consequences!\n")
        ##break() #out from for loop
      }


      crownBaseMissing <- TRUE
      crownBase <- NA


      ### READING CROWN BASE FROM FILE ####
      if(is.element("crownBase", colnames(nowMeta))){
        #TODO, never used the reading from file...
        crownBase <- nowMeta$crownBase
          if(!length(crownBase)==1){
            if(length(crownBase)==0){
              crownBase <- NA
            } else {
              crownBase <- crownBase[length(crownBase)]
              cat("Warning: more than one crown base value for tree",treeName,"- we take the last one.\n")
            }
          }
          cat("Input crown base for tree",treeName,"is at z",crownBase,"\n")

          crownStart <- crownBase-minZ
          crownLength <- treeLAS@header@PHB$`Max Z` - crownBase



          png(paste0(crownImagePath,treeName,"_crown_mes.png"), height = 800, width = 1400, type = "cairo")
          par(mfrow=c(2,2))
          if(do.Plot){
            plot(tempTree@data$Y, tempTree@data$Z, asp = 1, cex = 0.1)
            abline(h = c(minZ,minZ + crownStart), col = "red")
            abline(h = c(treeLAS@header@PHB$`Max Z`, minZ + crownStart+0.1), col = "black")
            #abline(h = (minZ + 1.3 - groundCutHeight), col = "green")

          }

          #metaVars$est.height <- round(crownLength+crownStart+groundCutHeight,2)
          metaVars$est.crownLength <- round(crownLength,2)
          metaVars$est.crownBase <- round(crownBase,2)


          if(do.Plot){
            plotTitle <- paste0(treeName, " - manual crown base at ", round(crownBase,1))
            if(is.element("dbh", colnames(nowMeta))){
              lines(c(nowMeta$y - 1, nowMeta$y + 1),c(nowMeta$z,nowMeta$z), col = "green")
              mtext(paste0("DBH=",round(nowMeta$dbh,1),"cm"), side = 2, adj = 0, col = "green")
              mtext(paste0("at z=",round(nowMeta$z,1),""), side = 2, adj = 1, col = "green")
            }
            mtext(paste0("stem=",round(crownStart,1),"m"), col = "red", side = 4, adj = 0)
            mtext(paste0("crown=",round(crownLength,1),"m (",round(crownLength/(crownLength+crownStart)*100),"%)"), side = 4, adj = 1)
            title(plotTitle)
          }

      } else
        {


        cb1 <- Sys.time()
        cat("Crown base detection ")
        # produces value of crownBase (has added 100 in it) and crownStart (relative to floor, no 100 in it)
        # and centers, a dataframe of center coordinates and diameters until crown start



        referenceCounter <- 0 #counts when to start taking reference diameter
        referenceDiameter <- 100 #initially 100m
        if(is.element("dbh", colnames(nowMeta))){
          try({
            referenceDiameter <- nowMeta$dbh/100
          })
        }
        refList <- vector()
        j <- 1

        png(paste0(stemAnalysisPath,"",treeName,"_stem_diameters.png"), height = 6000, width = 6000, type = "cairo")
        par(mfrow=c(ceiling(length(heightGrip)/20),20),oma = c(0, 0, 3, 0))

        ### CROWN BASE DETECTION LOOP ####
        previousCenter <- firstCenter

        dbhPosition <- which(heightGrip == 1.3)
        if(length(dbhPosition) == 0){
          dbhPosition <- 9999999
          measureOrder <- c(1:length(heightGrip))
        } else {
          measureOrder <- c(dbhPosition:1, (dbhPosition+1):length(heightGrip))
        }


        considerTrunkStart <- -1 # if slices down from dbh already
        # contain the crown, this value holds
        # how many slices we have to consider
        for(k in 1:length(measureOrder)){ #for every slice of stem diameter
          j <- measureOrder[k]
          # SLICING
          sliceLAS <- filter_poi(treeLAS, Z >= (minZ + heightGrip[j]), Z < (minZ + heightGrip[j] + slice.height))

          #sliceLAS <- filter_poi(sliceLAS, Classification != LASGROUND)

          if(k == dbhPosition + 1){
            # now going upwards
            previousCenter <- firstCenter # reset to DBH location
          }

          nextCircleRadiusIncrease <- 0.05 # 5 cm more diameter for next stem circle
          multiplyReferenceD <- 1
          if(j < 10){
            multiplyReferenceD <- (10-j)^2/(10^2)*2+1
          }
          nowClipRadius <- (previousCenter$d/2*multiplyReferenceD) + nextCircleRadiusIncrease



          if(sliceLAS@header@PHB$`Number of point records` == 0){
            if(do.Plot) plot(x = 0, y = 0, asp = 1, main = "",
                             type = "n", xlab = "x", ylab = "y")
          } else {

            if(crownBaseMissing || k <= dbhPosition){
              if(do.Plot) plot(sliceLAS@data$X, sliceLAS@data$Y,
                               asp = 1, main = "",
                               type = "n", xlab = "x", ylab = "y")
            } else {
              if(do.Plot) plot(0,
                               xlim = c(previousCenter$x - 3*nowClipRadius,
                                        previousCenter$x + 3*nowClipRadius),
                               ylim = c(previousCenter$y - 1.6*nowClipRadius,
                                        previousCenter$y + 1.6*nowClipRadius),
                               asp = 1, main = "",
                               type = "n", xlab = "x", ylab = "y")
            }


          }

          if(k <= dbhPosition){
            box(col="red")
            # going downwards to stem base indicated by a red box
          }

          if(sliceLAS@header@PHB$`Number of point records`<=20){
            if(sliceLAS@header@PHB$`Number of point records`==0){
              if(do.Plot) title(main = paste0("h=",heightGrip[j],"m (z=",round(minZ + heightGrip[j],2),") EMPTY"))
              #cat("Empty Slice at",(minZ + heightGrip[j]),"m - next!\n")
              next()
            }
            center <- data.frame("x" = median(sliceLAS@data$X, na.rm = T), "y" = median(sliceLAS@data$Y, na.rm = T),
                                 "z" = (minZ + heightGrip[j]), "d" = -1, "lower" = heightGrip[j])
            centers[j,] <- center
            if(do.Plot) points(sliceLAS@data$X, sliceLAS@data$Y, col = "green",  cex = 0.001)
            if(do.Plot) points(center$x, center$y, cex = 2, col = "red")
            if(do.Plot) title(main = paste0("h=",heightGrip[j],"m (z=",round(minZ + heightGrip[j],2),") LESS THAN 20 POINTS"))
            #cat("There are too few points in slice at",(minZ + heightGrip[j]),"m - next!\n")
            next()
          }




          stemLAS <- clip_circle(sliceLAS, previousCenter$x, previousCenter$y, nowClipRadius)




          if(crownBaseMissing){
            # FILTERING
            referenceCounter <- referenceCounter + 1
            if(do.Plot) points(sliceLAS@data$X, sliceLAS@data$Y,
                               col = "red", pch = ".")
            filter.x <- quantile(sliceLAS@data$X, c(filter.quant, 1 - filter.quant))
            filter.y <- quantile(sliceLAS@data$Y, c(filter.quant, 1 - filter.quant))
            sliceLAS <- filter_poi(sliceLAS, X > min(filter.x)-tolerance, X < max(filter.x)+tolerance,
                                   Y > min(filter.y)-tolerance, Y < max(filter.y)+tolerance)
            rm(filter.x, filter.y)
            if(do.Plot) points(sliceLAS@data$X, sliceLAS@data$Y, col = "blue", pch = ".")
            if(do.Plot) points(stemLAS@data$X, stemLAS@data$Y,
                               col = "black",  pch = ".")

            # FITTING DIAMETER CIRCLE
            # CircleFitByLandau used by now (tight grip)
            c2 <- CircleFitByLandau(cbind(sliceLAS@data$X, sliceLAS@data$Y))
            center <- data.frame("x" = c2[1], "y" = c2[2], "z" = (minZ + heightGrip[j]),
                                 "d" = c2[3]*2, "lower" = heightGrip[j])

            if(do.Plot) draw.circle(center$x, center$y, radius = center$d/2, border = "red", lty = 2)

            # criteria for crown basal height!
            exceedingLimitForCrownBase <- center$d > 2*referenceDiameter



            if(do.Plot) mtext(paste0("d.cr=",round(center$d*100,1),"cm  "),
                  side = 3, adj = 1, line = -2, col = "#2244CC")
            if(exceedingLimitForCrownBase){
              if(do.Plot) mtext("EXCEEDS  ",
                    side = 3, adj = 1, line = -3.5, col = "#2244CC")
            }

            centers[j,] <- center
            #print(center)



          } else {
            # plot crown points in red with right xlim values
            if(do.Plot) points(sliceLAS@data$X, sliceLAS@data$Y,
                             col = "red", pch = ".")
            # plot points of the slice in green
            if(do.Plot) points(stemLAS@data$X, stemLAS@data$Y,
                               col = "#448711", pch = ".")



          }

          if(do.Plot) mtext(paste0(" x=", round(previousCenter$x,2), "\n",
                                   " y=", round(previousCenter$y,2), ""),
                            side = 1, adj = 0, line = -1.5, col = "black")
          if(do.Plot) mtext(paste0("cl.d=",round(nowClipRadius*100,1)*2," \n",
                                   "prv.d=",round(previousCenter$d*100,1), " "),
                            side = 1, adj = 1, line = -1.5, col = "black")

          if(stemLAS@header@PHB$`Number of point records`> 10){
            #cat("Too few points in this stem circle at", (minZ + heightGrip[j]), "m - next!...\n")

            c3stem <- CircleFitByLandau(cbind(stemLAS@data$X, stemLAS@data$Y),
                                        ParIni = c(previousCenter$x, previousCenter$y, previousCenter$d/2))

            stemCenter <- data.frame("x" = c3stem[1], "y" = c3stem[2], "z" = (minZ + heightGrip[j]),
                                     "d" = c3stem[3]*2, "movedCM" = 0,
                                     "angle" = 0, "lower" = heightGrip[j])
            stemMovedAngle <- atan2((stemCenter$y - previousCenter$y),
                                    (stemCenter$x - previousCenter$x))/2/pi*360
            if(stemMovedAngle < -45){
              stemMovedAngle <- stemMovedAngle + 360
            }
            stemMovedDist <- sqrt((stemCenter$y - previousCenter$y)^2 +
                                    (stemCenter$x - previousCenter$x)^2)

            stemCenter$movedCM <- round(stemMovedDist*100,2)
            stemCenter$angle <- round(stemMovedAngle,1)



            flatNess <- (stemCenter$movedCM)/abs(stemCenter$z - previousCenter$z)

            if(do.Plot) mtext(paste0("  +",stemCenter$movedCM," cm\n",
                                     "  +", round(flatNess, 1),"%\n",
                                     "  ", stemCenter$angle,"°"),
                              side = 3, adj = 0, line = -4.5, col = "#448711")

            if(do.Plot) draw.circle(stemCenter$x, stemCenter$y,
                                    radius = stemCenter$d/2, border = "black", lty = 2, lwd = 2)

            if(do.Plot) title(main = paste0("h=",heightGrip[j],"m (z=",
                                            round(minZ + heightGrip[j],2),"), d=",
                                            round(stemCenter$d*100,1),"cm"))


            stems[j,] <- stemCenter

            if(stemCenter$d < (previousCenter$d - 0.1) |
               stemCenter$d > (previousCenter$d + 0.05) |
               flatNess > 150){

              #stemCenter$d <- previousCenter$d
              # if suddenly the diameter is shrinked by more than 10 cm,
              # or if suddenly expanded by more than 5 cm,
              # or moved more than 5 cm (k = 100%) each 5 cm slice
              # rather stay with the old diameter, seems more stable...
            } else {
              previousCenter <- stemCenter
            }



          } else {
            if(do.Plot) title(main = paste0("h=",heightGrip[j],"m (z=",
                                            round(minZ + heightGrip[j],2),"), LESS THAN 10 STEM POINTS"))

          }



          if(crownBaseMissing || k <= dbhPosition){
            # DETERMINING reference diameter for crown base
            # only if crown base was not found yet

            if(referenceCounter <= 10 && referenceCounter > 5){
              # reference Diameter is taken between 5x slice (25 cm) and 10x slice (50 cm) from floor
              # if slice of course is 5 cm...
              refList <- c(refList, center$d)
            }
            # we take referenceDiameter now from tree_dbh.txt list!
            # only if we don't have a valid referenceDiameter, then we take the newly measured one
            if(referenceDiameter == 100){
              # nowMeta doesn't contain a dbh or bhd field, so we measure it from slice 5 - 10
              if(referenceCounter == 10){
                referenceDiameter <- mean(refList)
                cat("- CALCULATING NEW REFERENCE STEM DIAMETER",round(referenceDiameter*100,2),"cm.\n")
                cat("CBH Crown basal height ")

                if(referenceDiameter > referenceDiameterLimit/100){
                  #exceededReferenceLimit <- exceededReferenceLimit + 1
                  #treeTooBig <- T
                  #cat("Too big reference diameter for",treeName," - we take the alternative ratio!\n")
                  cat("Very big reference diameter for",treeName," - no consequences!\n")
                }
              }
            }

            # CHECKING if we are now at crown base
            if(referenceCounter >= 0){ # change noiw 06.02.2025 crown base can be zero!
              # above 50 cms we can detect crown starting point (10*5cm)
              if(exceedingLimitForCrownBase){
                crownCounter <- crownCounter + 1
                if(crownCounter > 1/slice.height || #must be steady over 1 m if going down from dbh
                   (k > dbhPosition && # if going up from dbh, we need to consider the trunkstart below
                    (crownCounter + considerTrunkStart) > 1/slice.height)){
                  if(j < 1.3){
                    crownStart <- heightGrip[j] #go down again for that one meter (20 slices)
                  } else {
                    crownStart <- heightGrip[trunc(j - 1/slice.height)] #go down again for that one meter (20 slices)
                  }
                  crownBase <- crownStart+minZ
                  cat("at z =",crownBase,"(rel. height",crownStart,"m).\n")
                  crownBaseMissing <- FALSE

                  crownLength <- treeLAS@header@PHB$`Max Z` - crownBase

                  #cat("Automatized crown base for tree",treeName,"is at z",crownBase,"\n")
                  #break() #out from for loop
                  # new 2025-01-29 continue for stem circles upwards
                }
              } else {
                if(k <= dbhPosition){
                  if(considerTrunkStart == -1){
                    considerTrunkStart <- crownCounter
                  }
                } else {
                  if(considerTrunkStart > 0){
                    considerTrunkStart <- 0 #reset trunk start, if stem above dbh does not have crown points in it
                  }
                }
                # as soon as diameter springs back to normal, we need to reset crownbase counter
                crownCounter <- 0
              }
            }
          }
        } #for every slice of stem diameter
        #rm(referenceCounter, crownCounter, j, sliceLAS, center, refList, heightGrip, c2)

        mtext(paste0(treeName," stem and crown circles   reference Diameter: ",round(referenceDiameter*100,1),"cm    minZ = ",round(minZ,2),". "), outer = TRUE, cex = 2.0)


        # conclude image of slices, output and start new picture with crown base
        dev.off()



        centerPoint1m <- stems[trunc(1/slice.height),]
        centerPoint2m <- stems[trunc(2/slice.height),]
        centerPoint160cm <- stems[trunc(1.6/slice.height),]


        centers <- centers[!is.na(centers$d),]
        stems <- stems[!is.na(stems$d),]

        if(is.na(centerPoint1m$x)){
          centerPoint1m <- stems[trunc(1/slice.height),] # if first measurement was invalid
        }
        if(is.na(centerPoint2m$x)){
          centerPoint2m <- stems[trunc(2/slice.height),]
        }
        if(is.na(centerPoint160cm$x)){
          centerPoint160cm <- stems[trunc(2/slice.height),]
        }


        metaVars$center100cmX <- centerPoint1m$x
        metaVars$center100cmY <- centerPoint1m$y
        metaVars$center100cmZ <- centerPoint1m$z
        metaVars$center160cmX <- centerPoint160cm$x
        metaVars$center160cmY <- centerPoint160cm$y
        metaVars$center160cmZ <- centerPoint160cm$z
        metaVars$center200cmX <- centerPoint2m$x
        metaVars$center200cmY <- centerPoint2m$y
        metaVars$center200cmZ <- centerPoint2m$z
        metaVars$center100cmD <- centerPoint1m$d
        metaVars$center160cmD <- centerPoint160cm$d
        metaVars$center200cmD <- centerPoint2m$d

        metaVars$stemLeanDist1to2 <- sqrt((centerPoint2m$y - centerPoint1m$y)^2 + (centerPoint2m$x - centerPoint1m$x)^2)
        metaVars$stemLeanAngle1to2 <- atan2((centerPoint2m$y - centerPoint1m$y), (centerPoint2m$x - centerPoint1m$x))/2/pi*360

        metaVars$stemLeanDist1to160 <- sqrt((centerPoint160cm$y - centerPoint1m$y)^2 + (centerPoint160cm$x - centerPoint1m$x)^2)
        metaVars$stemLeanAngle1to160 <- atan2((centerPoint160cm$y - centerPoint1m$y), (centerPoint160cm$x - centerPoint1m$x))/2/pi*360

        metaVars$taper1to2 <- centerPoint1m$d - centerPoint2m$d
        metaVars$taper1to160 <- centerPoint1m$d - centerPoint160cm$d

        write.table(centers, file = paste0(stemAnalysisPath,treeName,"_crown_centers.txt"), row.names = FALSE, sep = "\t", dec = ".")
        write.table(stems, file = paste0(stemAnalysisPath,treeName,"_stem_centers.txt"), row.names = FALSE, sep = "\t", dec = ".")

        if(crownBaseMissing){
          ## NO CROWN DETECTED SECTION ###
          crownBase <- altCrownBase
          crownStart <- crownBase-minZ
          cat("not detected! Estimate", alternativeCrownBase.Ratio*100, "% to z=", crownBase,"m.\n")
          crownLength <- treeLAS@header@PHB$`Max Z` - crownBase

          metaVars$est.crownAreaAlt <- -1
          #metaVars$est.height <- round(crownLength+crownStart+groundCutHeight,2)
          metaVars$est.crownLength <- round(crownLength,2)
          metaVars$est.crownBase <- round(crownBase,2)

        }


        {

          png(paste0(crownImagePath,treeName,"_crown.png"), height = 800, width = 1400, type = "cairo")
          par(mfrow=c(2,2))
          if(do.Plot){
            plot(tempTree@data$Y, tempTree@data$Z, asp = 1, cex = 0.1)
            segments(x0 = centers[trunc(2/slice.height),]$y, y0 = minZ, #2/slice.height refers to 2m according to slices
                     x1 = centers[trunc(j-1/slice.height)-1,]$y, y1 = minZ + crownStart, col = "red")
            abline(h = c(minZ,minZ + crownStart), col = "red")
            abline(h = c(treeLAS@header@PHB$`Max Z`, minZ + crownStart+0.1), col = "black")
            #abline(h = (minZ + 1.3 - groundCutHeight), col = "green")

          }

          DBH <- centers[centers$lower==1.3-groundCutHeight,]
          # DBH <- DBH[!is.na(DBH$d),]
          # metaVars$est.DBH <- round(DBH$d*100,1)
          # metaVars$est.x.DBH <- DBH$x
          # metaVars$est.y.DBH <- DBH$y
          # metaVars$est.z.DBH <- nowMeta$z
          # cat(paste0("DBH alternative ",round(DBH$d*100,1)," cm at position x="
          #            ,round(DBH$x,2),"  y=",round(DBH$y,2),"  z=",round(nowMeta$z,2),"m\n"))
          cat("HERE DBH ALTERNATIVE IS NOT IMPLEMENTED YET!\n")

          #metaVars$est.height <- round(crownLength+crownStart+groundCutHeight,2)
          metaVars$est.crownLength <- round(crownLength,2)
          metaVars$est.crownBase <- round(crownBase,2)

          if(do.Plot){
            if(crownBaseMissing){
              plotTitle <- paste0(treeName, " - cbh at ", alternativeCrownBase.Ratio*100 ,"%")
            } else {
              plotTitle <- treeName
            }
            if(is.element("dbh", colnames(nowMeta))){
              lines(c(nowMeta$y - 1, nowMeta$y + 1),c(nowMeta$z,nowMeta$z), col = "green")
              mtext(paste0("DBH=",round(nowMeta$dbh,1),"cm"), side = 2, adj = 0, col = "green")
              mtext(paste0("at z=",round(nowMeta$z,1),""), side = 2, adj = 1, col = "green")
            }
            mtext(paste0("stem=",round(crownStart,1),"m"), col = "red", side = 4, adj = 0)
            mtext(paste0("crown=",round(crownLength,1),"m (",round(crownLength/(crownLength+crownStart)*100),"%)"), side = 4, adj = 1)
            title(plotTitle)


            #TODO TREELIST ADD DBH AND HEIGHT COLUMN IF EXISTS FROM INVENTORY
            if(FALSE){
              # compare with DBH and height from inventory list
              DBH_bias <- round((DBH$d*100 - infoTreesMatched$DBH[i])/infoTreesMatched$DBH[i]*100,1)
              crownBase_bias <- round((crownBase - infoTreesMatched$crownBase[i])/infoTreesMatched$crownBase[i]*100,1)
              abline(h = infoTreesMatched$crownBase[i]+minZ, col = "blue", lty = "dashed")

              h_auto <- treeLAS@header@PHB$`Max Z` - minZ + groundCutHeight
              height_bias <- round((h_auto - infoTreesMatched$height[i])/infoTreesMatched$height[i]*100,1)



              legend("bottomleft", c(paste0(infoTreesMatched$species[i]," = ",treeSpecies(infoTreesMatched$species[i], TRUE)),
                                     #legend("bottomleft", c(paste0(numberSp," = ",species),
                                     paste0("DBH=",round(DBH$d*100,1),"cm (ref=",round(infoTreesMatched$DBH[i],1),
                                            "cm, bias=",DBH_bias,"%)"),
                                     paste0("height=",round(h_auto,1),"m (ref=",round(infoTreesMatched$height[i],1),
                                            "m, bias=",height_bias,"%)"),
                                     paste0("crb-height=",round(crownBase,1),"m (ref=",round(infoTreesMatched$crownBase[i],1),
                                            "m, bias=",crownBase_bias,"%)")),
                     text.col = c("black","green","red","blue"), bty="n", cex = 1.6, inset = c(-.03, .02))



              mtext(paste0("stem=",round(crownStart,1),"m"), col = "red", side = 4, adj = 0)
              mtext(paste0("crown=",round(crownLength,1),"m (",round(crownLength/(crownLength+crownStart)*100),"%)"), side = 4, adj = 1)


            }

          }

          cb2 <- Sys.time()
          timeNB <- as.difftime(cb2 - cb1)
          #cat("End of crown base searching - in", round(timeNB,1), units(timeNB),"\n")
        }







      }





      # HULLING crown base and getting volume ####
      # generates a alphahull3d over the crown (crownHull3D)
      #           a alphahull over crown xy view (crownOutlineXY)

      cat("CV Crown hull volume ")
      {
        start <- Sys.time()

        #metaVars$est.crownVolume <- 1

        crownLAS <- filter_poi(treeLAS, Z >= crownBase)
        distVect <- sqrt((crownLAS@data$X - nowMeta$x)^2 + (crownLAS@data$Y - nowMeta$x)^2)

        metaVars$meanCrownX <- mean(crownLAS@data$X)
        metaVars$meanCrownY <- mean(crownLAS@data$Y)
        metaVars$meanCrownZ <- mean(crownLAS@data$Z)
        metaVars$meanweightedCrownX <- weighted.mean(crownLAS@data$X, distVect)
        metaVars$meanweightedCrownY <- weighted.mean(crownLAS@data$Y, distVect)
        metaVars$meanweightedCrownZ <- weighted.mean(crownLAS@data$Z, distVect)


        # all 15x15 cm we keep 3 points
        set.seed(15)
        if(detail.level < 0){
          crownLAS <- decimate_points(crownLAS, random_per_voxel(res = round(-detail.level,1), n = 5))
        } else if(detail.level == 0){
          crownLAS <- decimate_points(crownLAS, random_per_voxel(res = 0.25, n = 2))
        } else if(detail.level > 0){
          crownLAS <- decimate_points(crownLAS, random_per_voxel(res = 0.10, n = round(detail.level)))
        }


        if(crownLAS@header@PHB$`Number of point records` > 300000){
          #exceededReferenceLimit <- exceededReferenceLimit + 1
          cat("More than 300000 crown points - it is a big tree (no consequences)!\n")
          ##break() #out from for loop
        }


        #metaVars$est.crownVolume <- 2
        #plot(crownLAS)
        matrix <- cbind("X" = crownLAS@data$X,
                        "Y" = crownLAS@data$Y,
                        "Z" = crownLAS@data$Z)

        #metaVars$est.crownVolume <- trunc(crownLAS@data$X[1])
        crownHull3D <- ashape3d(matrix, alpha = vol.alpha)
        # volume_ashape3d(crownHull3D)
        metaVars$est.crownVolume <- round(volume_ashape3d(crownHull3D),
                                          1)
        cat("is", metaVars$est.crownVolume, "m3 - plot ")
        #plot(crownHull3D)

        triangs <- crownHull3D$triang[crownHull3D$triang[,9]==2,] #got this line by accident, says only outer lines of net

        #metaVars$est.crownVolume <- length(triangs)

        mainNow <- treeName

        mainNow <- paste0(treeName, " - ", detail.level, " detail")
        # if(bigTree){
        #   mainNow <- paste0(treeName, " - BIG")
        # }
        plot(crownLAS@data$Y, crownLAS@data$Z, asp = 1, cex = 0.5, main = mainNow)


        shapeHullPoints <- data.frame()

        x_axis <- 2 # 1 x  2 y   3 z
        y_axis <- 3 # 1 x  2 y   3 z
        for(a in 1:length(triangs[,1])){
          p1 <- matrix[triangs[a,1],]
          p2 <- matrix[triangs[a,2],]
          p3 <- matrix[triangs[a,3],]
          newPoints <- data.frame(rbind(p1[1:2], p2[1:2], p3[1:2]))
          shapeHullPoints <- rbind(shapeHullPoints, newPoints)
          segments(x0 = p1[x_axis], y0 = p1[y_axis], x1 = p2[x_axis], y1 = p2[y_axis], col = "red")
          segments(x0 = p1[x_axis], y0 = p1[y_axis], x1 = p3[x_axis], y1 = p3[y_axis], col = "red")
          segments(x0 = p3[x_axis], y0 = p3[y_axis], x1 = p2[x_axis], y1 = p2[y_axis], col = "red")
        }
        rm(x_axis, y_axis, p1, p2, p3, newPoints, matrix)
        rm(triangs)
        gc()
        mtext(paste0("Crown volume: \n",round(volume_ashape3d(crownHull3D),1)," m3."), side = 4, adj = 1)


        #polygon(borderPoints$x, borderPoints$y)
        #colnames(shapeHullPoints) <- c("x", "y")
        shapeHullPoints <- (unique(shapeHullPoints))
        #crownOutlineXY <- ahull(shapeHullPoints$x, shapeHullPoints$y, alpha = area.alpha)

        # plot(tempTree@data$X, tempTree@data$Y, asp = 1, cex = 0.1, main = treeName)

        stop <- Sys.time()
        timeNB <- as.difftime(stop - start)
        cat("done in",round(timeNB,1),units(timeNB),"\n")




        cat("CA Crown projection area ")
        {
          plot.now <- TRUE
          t1 <- Sys.time()

          if(do.Plot) tempTree <- crownLAS

          crownLAS@data$Z <- 0
          crownLAS <- decimate_points(crownLAS, homogenize(density = 30, res = 0.2))
          crownLAS <- filter_duplicates(crownLAS)

          crownOutlineXY <- ahull(crownLAS@data$X, crownLAS@data$Y, alpha = area.alpha)
          metaVars$est.crownArea <- round(areaahulleval(crownOutlineXY),1)
          cat("is", metaVars$est.crownArea, "m2 - plot ")
          t2 <- Sys.time()
          #print.difftime(t2-t1)

          if(do.Plot)  plot(tempTree@data$Y ~ tempTree@data$X, cex = 0.0001, asp = 1)
          if(do.Plot)  points(crownLAS@data$Y ~ crownLAS@data$X, cex = 0.5, pch = "+", asp = 1, col = "red")
          if(do.Plot)  plot(crownOutlineXY, add = T, col = "red", wpoints = F)
          if(do.Plot) mtext(paste0("Crown area: \n",metaVars$est.crownArea," m2."),
                            side = 4, adj = 1)

        }




        timeNB <- as.difftime(Sys.time() - start)
        cat("done in",round(timeNB,1),units(timeNB),"\n")
        #cat("Crown hulling in",round(timeNB,1),units(timeNB),"\n")
      }


      # Drawing image of diagonals (min max width of furthest point) ####
      # generates two diagonal info fields (diag.max - length and diag.max.angle also for min)
      #           and a list of outer crown projection points (borderPoints)
      if(1==2){
        start <- Sys.time()
        borderPoints <- xyset[ashape(xyset$x, xyset$y, alpha = 0.3)$alpha.extremes,]
        if(do.Plot) plot(borderPoints, asp = 1, col = "blue", main = paste0("Diagonals of ",treeName))
        md <- as.matrix(dist(borderPoints))
        maxs <- apply(md, 1, max)

        linePoints <- which(md == max(maxs), arr.ind = TRUE)[1,] #maximum diagonal
        diag.max <- max(maxs)
        diag.max.angle <- round(atan((borderPoints[linePoints[1],2]-borderPoints[linePoints[2],2])/
                                       (borderPoints[linePoints[1],1]-borderPoints[linePoints[2],1]))*180/pi,2)+90
        borderPoints[linePoints[1],]
        borderPoints[linePoints[2],]

        if(do.Plot) {
          segments(x0 = borderPoints[linePoints[1],1], y0 = borderPoints[linePoints[1],2],
                   x1 = borderPoints[linePoints[2],1], y1 = borderPoints[linePoints[2],2], col = "red")
          mtext(paste0("max=",round(diag.max,1),"m ",round(diag.max.angle),"deg"), side = 3, adj = 0, col = "red")
        }


        linePoints <- which(md == min(maxs), arr.ind = TRUE)[1,] #minimum diagonal of furthest points
        diag.min <- min(maxs)
        diag.min.angle <- round(atan((borderPoints[linePoints[1],2]-borderPoints[linePoints[2],2])/
                                       (borderPoints[linePoints[1],1]-borderPoints[linePoints[2],1]))*180/pi,2)+90
        borderPoints[linePoints[1],]
        borderPoints[linePoints[2],]


        if(do.Plot) {
          segments(x0 = borderPoints[linePoints[1],1], y0 = borderPoints[linePoints[1],2],
                   x1 = borderPoints[linePoints[2],1], y1 = borderPoints[linePoints[2],2], col = "orange")
          mtext(paste0("min=",round(diag.min,1),"m ",round(diag.min.angle),"deg"), side = 3, adj = 1, col = "orange")


          mtext(paste0("Crown area: \n",round(areaahulleval(crownOutlineXY),1)," m2.\n"), side = 1, adj = 1)
          mtext(paste0("Crown volume: \n",round(volume_ashape3d(crownHull3D),1)," m3.\n"), side = 1, adj = 0)

          #points(centers[20:50,1:2], col = "green", pch = 4)
        }

        # diagonals meaning: starting from north turning counterclockwise (positive)
        #  0 is N-S     45 is NW-SE     90 is W-E    135 is SW-NE and 178 is SSW-NNE

        crownCircle <- CircleFitByLandau(borderPoints)
        draw.circle(x = crownCircle[1], y = crownCircle[2], radius = crownCircle[3], border = "lightblue")
        points(x = crownCircle[1], y = crownCircle[2], col = "lightblue", pch = 4, cex = 2)





        metaVars$est.crownVolume <- round(volume_ashape3d(crownHull3D),1)

        metaVars$est.DiamCrownMax <- round(diag.max,1)
        metaVars$est.DiamCrownMax.Angle <- round(diag.max.angle,0)
        metaVars$est.DiamCrownMin <- round(diag.min,1)
        metaVars$est.DiamCrownMin.Angle <- round(diag.min.angle,0)

        metaVars$est.crownDiameter <- round(crownCircle[3]*2,1)
        metaVars$est.x.crownDiameter <- crownCircle[1]
        metaVars$est.y.crownDiameter <- crownCircle[2]





        rm(maxs, linePoints, xyset)
        timeNB <- as.difftime(Sys.time() - start)
        cat("Diagonals drawing in",round(timeNB,1),units(timeNB),"\n")
      }
    }



    gstop <- Sys.time()
    timeNB <- as.difftime(gstop - gstart)
    mtext(paste0("Time ",round(timeNB,1),units(timeNB),"\n\n\n\n"), side = 1, adj = 1)
    dev.off()
    cat("This single tree took:",round(timeNB,1),units(timeNB),"\n\n")


  }, error = function(e) {
    # If Error happens at one tree
    print(paste("MY_ERROR:  ",e))
    cat("Error catching - closing jpeg device...\n\n")
    try(dev.off())

  })


  return(metaVars)
}


#' Calculation of individual tree measurements
#'
#' Input files are  single tree .laz files, there is no splitting of
#' the segmented crowns file.
#'
#' @export
computeCrownParams <- function(fileFinder, loopStart = 1, loopEnd = 0,
                               drawImages = T,
                               remeasureDBH = F,
                               decimateTreeForFasterCrowns = FALSE,

                               mode = "ALLGO", maxRadius = 0, ipad = FALSE,

                               cutWindow = c(-1000,-1000,2000),
                               zScale = 2, limitShare = 0.004, voxelSize = 0,

                               vol.alpha = 2, alternativeCrownBase.Ratio = 0.3,
                               fogFilter.estHeight = FALSE,
                               detail.level = 0,

                               selector = "xyzcit0",
                               limitSpanSide = 30, limitSpanArea = 440,
                               referenceDiameterLimit = 150,

                               quantileIntensity = 15, numberIntensity = 0, CC_level = 10, CC_numberPoints = 1000,
                               clipHeight = 3, bottomCut = 1, bushPreparation = FALSE, filterSOR = FALSE,

                               do.CrownBase = TRUE,
                               groundCutHeight = 0, silent = TRUE, retainPointClouds = FALSE, nr_cores = 0,


                               treeNames = NULL,
                               smallTrees = FALSE,
                               dirPath = paste0(getwd(),"/")){

  {
    if(ipad){
      clipHeight <- 2.2
      bottomCut <- 0.2
      selector <- "xyzcitRGB0"
    }
    gstart <- Sys.time()

    fileFinder <- removeUmlaut(fileFinder)



    if(!exists("outLAS")){
      outLAS <<- NA # all crowns together in one las file after crownFeeling()
      outLAS_name <<- ""
    }


    if(mode == "COMP"){
      setString <- generateSetString(fileFinder, mode = mode,
                                     threshold = numberIntensity, threshold_percent = quantileIntensity,
                                     bottomCut = bottomCut, clipHeight = clipHeight,
                                     bushPreparation = bushPreparation, filterSOR = filterSOR,
                                     level = CC_level, numberOfPoints = CC_numberPoints, silent = silent)
      dbhPath <- paste0(dirPath,setString,"_dbh/")

      if(!dir.exists(dbhPath)){
        setString <- generateSetString(fileFinder, mode = mode, cutWindow = cutWindow,
                                       threshold = numberIntensity, threshold_percent = quantileIntensity,
                                       bottomCut = bottomCut, clipHeight = clipHeight,
                                       bushPreparation = bushPreparation, filterSOR = filterSOR,
                                       level = CC_level, numberOfPoints = CC_numberPoints, silent = silent)
        dbhPath <- paste0(dirPath,setString,"_dbh/")
      }

    } else
      if(mode == "ALLGO") {
        dbhPath <- generateSetString(fileFinder = fileFinder, mode = mode,
                                     clipHeight = clipHeight, bottomCut = bottomCut,
                                     bushPreparation = bushPreparation,
                                     filterSOR = filterSOR, silent = TRUE)
        dbhPath <- paste0(dirPath,dbhPath,"_dbh/")

        if(!dir.exists(dbhPath)){
          dbhPath <- generateSetString(fileFinder = fileFinder, mode = mode, cutWindow = cutWindow,
                                       clipHeight = clipHeight, bottomCut = bottomCut,
                                       bushPreparation = bushPreparation,
                                       filterSOR = filterSOR, silent = TRUE)
          dbhPath <- paste0(dirPath,dbhPath,"_dbh/")
        }

        if(!dir.exists(dbhPath)){
          warning("There are no clustered input files! \nPlease run functions extractVegetation and clusterSplit/stemSplit before!\n")
          cat("Terminating on",format(Sys.time()),"\n\n")
        }


      }

    crownString <- generateSetString(limitShare = limitShare, zScale = zScale, voxelSize = voxelSize)
    crownPath <- paste0(dbhPath,"crowns",crownString,"/")
    setString_out <- crownString


    locationStr <- generateSetString(cutWindow = cutWindow, silent = silent)
    fileName <- paste0(fileFinder,"_crowns", locationStr,".laz")

    crownPathSingles <- paste0(crownPath,"singleTrees/")
    if(!dir.exists(crownPathSingles)) dir.create(crownPathSingles)


    #crownPathSingles <- paste0(crownPath,"singleTrees",format(Sys.time(), "%Y%m%d_%H%M"),"/") # nice idea timestamp on folder
    crownPathSinglesRef <- paste0(crownPath,"compareQuality/") #not used in this function







  }

  sink(paste0(crownPath,fileFinder,"_computeCrownParams_",format(Sys.time(), "%Y%m%d_%H%M"),"_Rcons.txt"), append = TRUE, split = TRUE)
  cat("COMPUTATION of Crown Parameters for case",fileFinder,"\n")
  cat("Job started:",format(gstart, "%Y-%m-%d, %H:%M"),"\n\n")


  cat("Collecting information from tree segmentation:\n")
  cat("   trees_dbh.txt... ")
  stemCoord <- read.table(paste0(dbhPath,"trees_dbh.txt"), header = TRUE)
  stemCoord$dist <- sqrt(stemCoord$x^2 + stemCoord$y^2)
  n.dbhTrees <- nrow(stemCoord)
  cat(n.dbhTrees, "trees found.\n")
  #cat("   trees_height.txt... ")
  #metaList <- read.table(paste0(dbhPath, "/trees_height.txt"), header = T)
  #n.hTrees <- nrow(metaList)
  #cat(n.hTrees, "trees found.\n")
  # loopEnd <- n.hTrees

  segmented.trees <- list.files(crownPathSingles, pattern = ".laz")
  cat("   laz in folder... ")
  cat(length(segmented.trees), "trees found.\n")


  if(loopEnd == 0 || loopEnd > n.dbhTrees){
    loopEnd <- n.dbhTrees
  }

  if(loopEnd < n.dbhTrees || loopStart != 1){
    cat("\n")
    cat("Subset all trees from", loopStart, "to", loopEnd, "\n")
    cat("The trees to measure are:\n   ")
    cat(stemCoord$id[c(loopStart:loopEnd)], sep = ", ")
    cat("\n\n")
  }

  ## also read here the crownBaseList!
  if(!do.CrownBase){
    crownBaseList.file <- paste0(crownPath, "/crownBaseList.txt")
    if(file.exists(crownBaseList.file)){
      crownBaseList <- read.table(crownBaseList.file, header = T, sep = "\t")
      colnames(crownBaseList)[1] <- "id"
      colnames(crownBaseList)[2] <- "crownBase"
    } else {
      cat("\nWARNING: There was no crownBaseList found!\n")
      do.CrownBase <- TRUE
    }
    #cat("USING OLD CROWN BASE FILE NOT IMPLEMENTED YET!\n")
  }


  #cat("Now go parallel and compute crown measures!\n")



  #cat("For now we go serial and calculate all tree parameters.\n\n")


  if(nr_cores == 0){
    cat("Automatic number of core setting:\n   ")
    nr_cores <- detectCores()
    cat(nr_cores, "cores available - take two less ")
    nr_cores <- nr_cores - 2
    if(nr_cores > 15){
      cat("but maximum 15!")
      nr_cores <- 15
    }

    # maxCores <- floor(totalRAM / oSize / 2) - 1
    # cat(maxCores, "cores recommended (RAM ratio of", round(totalRAM / oSize, 1), "/ 2 - 1).\n")
    # if(maxCores < 1){
    #   maxCores <- 1
    # }
    # if(maxCores <=  nr_cores){
    #   nr_cores <- maxCores
    # }

    cat("\n")
  }




  if(nr_cores == 1){
    # SERIAL SCRIPT

    cat("Serial measuring of", (loopEnd-loopStart+1), "trees:\n\n")
    timeStarted <- Sys.time()
    metaFrame <<- data.frame()

    for(i in c(loopStart:loopEnd)){

      if(is.null(treeNames)){
        #treeName <- paste0("tree",sprintf("%04d", trees[i]),".laz")
        treeName <- paste0(stemCoord$id[i])
        tryCatch({
          treeName <- paste0(sprintf("%04d", stemCoord$id[i]),"")
        }, error = function(e) {
          #cat("Warning - there are weird names present...\n")
        })
        if(treeName == "") treeName <- paste0(stemCoord$id[i])

      } else {
        treeName <- treeNames[i]
      }

      cat(i,"/",loopEnd,"- Starting with tree",treeName,
          "on",format(Sys.time(), "%Y-%m-%d, %H:%M"),"\n")

      nowMeta <- stemCoord[i,]

      if(!do.CrownBase){
        nowID <- stemCoord$id
        if(is.element(nowID, crownBaseList$id)){
          nowMeta$crownBase <- crownBaseList$crownBase[crownBaseList$id == nowID]
        }
      }


      # if(nowMeta$dist >= maxRadius & maxRadius != 0){
      #   writeLAZ = F
      #   writePicture = F
      # }

      outPut <- data.frame()

      nowTreeFile <- paste0(crownPathSingles, treeName, ".laz")

      if(file.exists(paste0(nowTreeFile))){
        tryCatch({
          outPut <- computeTree_i(treeLAS.path = nowTreeFile,
                                  decimateTreeForFasterCrowns = decimateTreeForFasterCrowns,
                                  treeName = treeName,
                                  fileFinder = fileFinder,
                                  drawImage = drawImages,
                                  detail.level = detail.level,
                                  nowMeta = nowMeta,
                                  remeasureDBH = remeasureDBH,
                                  measurePath = crownPath)
        })

        metaFrame <<- rbind(metaFrame, outPut)
      }


    }


    cat("Serial work done in a ")
    timeEnded <- Sys.time()
    print.difftime(round(timeEnded - timeStarted, 1))
    cat("\n")


  } else {
    # PARALLEL SCRIPT


    cl <- makeCluster(nr_cores, outfile="")
    registerDoParallel(cl)

    cat("Parallel measuring now of", (loopEnd-loopStart+1), "trees on", nr_cores,"cores...\n")
    i <- 1
    timePar1 <- Sys.time()

    fdc <<- foreach(i = loopStart:loopEnd, .errorhandling = 'remove',
                    .packages = c("treeX", "mgcv"), .verbose = FALSE) %dopar% {

                      # cat("Serial measuring of", (loopEnd-loopStart+1), "trees:\n\n")

                      # timeStarted <- Sys.time()
                      # metaFrame <<- data.frame()
                      # for(i in c(loopStart:loopEnd)){

                      if(is.null(treeNames)){
                        #treeName <- paste0("tree",sprintf("%04d", trees[i]),".laz")
                        treeName <- paste0(stemCoord$id[i])
                        tryCatch({
                          treeName <- paste0(sprintf("%04d", stemCoord$id[i]),"")
                        }, error = function(e) {
                          #cat("Warning - there are weird names present...\n")
                        })
                        if(treeName == "") treeName <- paste0(stemCoord$id[i])

                      } else {
                        treeName <- treeNames[i]
                      }

                      cat(i,"/",loopEnd,"- Starting with tree",treeName,
                          "on",format(Sys.time(), "%Y-%m-%d, %H:%M"),"\n")

                      nowMeta <- stemCoord[i,]


                      # if(nowMeta$dist >= maxRadius & maxRadius != 0){
                      #   writeLAZ = F
                      #   writePicture = F
                      # }

                      outPut <- data.frame()

                      nowTreeFile <- paste0(crownPathSingles, treeName, ".laz")

                      if(file.exists(paste0(nowTreeFile))){
                        tryCatch({
                          outPut <- computeTree_i(treeLAS.path = nowTreeFile,
                                                  decimateTreeForFasterCrowns = decimateTreeForFasterCrowns,
                                                  treeName = treeName,
                                                  fileFinder = fileFinder,
                                                  drawImage = drawImages,
                                                  detail.level = detail.level,
                                                  remeasureDBH = remeasureDBH,
                                                  nowMeta = nowMeta,
                                                  measurePath = crownPath)

                        })

                      }


                      #metaFrame <<- rbind(metaFrame, outPut)

                      return(outPut)
                    }


    #cat("Serial work done in a ")
    #timeEnded <- Sys.time()
    #print.difftime(timeEnded - timeStarted)
    stopCluster(cl)
    getDoParWorkers()

    cat("Parallel work done in a ")
    timePar2 <- Sys.time()
    print.difftime(round(timePar2 - timePar1, 1))
    cat("\n")

    # if(exceededReferenceLimit != 0){
    #   warning(paste0("MULTI-STEM: Skipping handling of", exceededReferenceLimit, "trees,",
    #                  "exceeding ref.diameter ", referenceDiameterLimit, " cm or span ", limitSpanSide, " m or area ", limitSpanArea, "m2.\n",
    #                  "Supposing border tree with unclear segmentation.\n",
    #                  "Assume crown base height at", alternativeCrownBase.Ratio*100, "% of tree height.\n",
    #                  "Set forceBigTrees true if you miss some trees!"))
    # }



    skippedTrees <- sum(lengths(fdc) == 0)
    if(skippedTrees > 0){
      cat("There were", skippedTrees, "trees out of the circle area.\n")

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
      for(j in c(3:length(df[1,]))){ # all num unless file and treeName
        df[,j] <- as.numeric(as.character(df[,j]))
      }
      #special case, convert treeName if it works
      tempTreeNames <- as.numeric(as.character(df[,2]))
      if(sum(is.na(tempTreeNames))==0){
        df[,2] <- tempTreeNames
      }
    }

    metaFrame <<- df



  }


  metaFrame$id <- (metaFrame$treeName)
  tryCatch({
    tempIDs <- strtoi(metaFrame$treeName, base = 10L)
    if(sum(is.na(tempIDs))<length(tempIDs)){
      metaFrame$id <- tempIDs
    } else {
      metaFrame$id <- metaFrame$treeName
    }
  }, error = function(e) {
    cat("Error in converting tree names to numbers - BEWARE!\n")
  })


  checkList <- (merge(metaFrame, stemCoord, by = "id", suffixes = c("", ".old")))
  write.table(checkList, file = paste0(crownPath, fileFinder, "_metaCrowns_all_",format(Sys.time(), "%Y%m%d_%H%M"),".txt"),
              sep = "\t", na = "", row.names = FALSE)


  #tempTreeNames <- tryCatch(as.integer(checkList$id), error=function(e) checkList$id, warning=function(w) checkList$id)
  tempTreeNames <- strtoi(checkList$id, base = 10L)
  if(sum(is.na(tempTreeNames)) > 0){
    tempTreeNames <- checkList$id
    }
  checkList$id <- tempTreeNames

# already done above
  if(remeasureDBH){
    applyThose <- !(checkList$est.x.DBH == 0 & checkList$est.y.DBH ==0)
    checkList$dbh[applyThose] <- checkList[applyThose,which(colnames(checkList) == "est.DBH")]
    checkList$x[applyThose]   <- checkList[applyThose,which(colnames(checkList) == "est.x.DBH")]
    checkList$y[applyThose]   <- checkList[applyThose,which(colnames(checkList) == "est.y.DBH")]
  }
  # long output format for crown parameters also measured
  checkList <- checkList[,c(which(colnames(checkList) == "x"),
                            which(colnames(checkList) == "y"),
                            which(colnames(checkList) == "z"),
                            which(colnames(checkList) == "id"),
                            which(colnames(checkList) == "name"),
                            which(colnames(checkList) == "cluster"),
                            which(colnames(checkList) == "randomCol"),
                            which(colnames(checkList) == "nPoints"),
                            which(colnames(checkList) == "dbh"),
                            which(colnames(checkList) == "d.circ"),
                            which(colnames(checkList) == "d.circ2"),
                            which(colnames(checkList) == "d.ell"),
                            which(colnames(checkList) == "d.gam"),
                            which(colnames(checkList) == "tegam.d_gam"),
                            which(colnames(checkList) == "est.height"),
                            which(colnames(checkList) == "est.crownBase"),
                            which(colnames(checkList) == "est.crownLength"),
                            which(colnames(checkList) == "est.crownArea"),
                            which(colnames(checkList) == "est.crownAreaAlt"),
                            which(colnames(checkList) == "est.crownVolume"),
                            #which(colnames(checkList) == "est.crownDiameter"),
                            #which(colnames(checkList) == "ld100circ"),
                            #which(colnames(checkList) == "ld200circ"),
                            which(colnames(checkList) == "species"),
                            which(colnames(checkList) == "inside"),
                            which(colnames(checkList) == "time_secs"))
  ]

  write.table(checkList, file = paste0(dbhPath, "trees_measured.txt"),
              sep = "\t", na = "", row.names = FALSE)




  gstop <- Sys.time()
  cat("Job ended:",format(Sys.time(), "%Y-%m-%d, %H:%M"),"\n\n")
  print.difftime(round(gstop - gstart, 1))
  gc()
  sink()


}

#' Calculation of individual tree measurements
#'
#' Loads in all single seperated Tree files from segmented _crowns.laz
#' and calculates for each individual the crown specific parameters.
#' E.g.: tree height, crown base, crown volume, crown projection area, ...
#' makes .png pictures and saves the list as
#' trees_height.txt and trees_measured.txt (incl. crown parameters)
#'
#' @param vol.alpha defines the alpha hull of crown volume, if 2 very rough (balloons) if 0.4 very tight (old settings)
#'
#' @param fileName origin of the segmented crown las file in dirPath
#' @param writeLAZ export the single trees as .laz
#' @param writePicture export the single trees as front, side and top view image, INCREASES COMPUTATION TIME
#' @param groundCutHeight how much distance between lowest point of tree and ground (if ground was cut off)
#' @param fogFilter.estHeight if TRUE the height will be calculated to 99.9 % quantile to the highest point, if FALSE to max height of the tree
#' @param maxRadius for Ebensee implemented, in a circle of 20 m all trees will be exported, but not the outer trees
#'
#' @export
computeTreeParams <- function(fileFinder, loopStart = 1, loopEnd = 0, getRAM = FALSE,
                              detail.level = 0,

                              drawNumberPoints = TRUE,

                              writeLAZ = TRUE, writePicture = TRUE, crownParameters = TRUE,

                              mode = "ALLGO", maxRadius = 0, ipad = FALSE,
                              cutWindow = c(-1000,-1000,2000), zScale = 2, limitShare = 0.004, voxelSize = 0,
                              vol.alpha = 2, alternativeCrownBase.Ratio = 0.3, fogFilter.estHeight = FALSE,

                              selector = "xyzcit0",
                              limitSpanSide = 30, limitSpanArea = 440,
                              referenceDiameterLimit = 150,

                              quantileIntensity = 15, numberIntensity = 0, CC_level = 10, CC_numberPoints = 1000,
                              clipHeight = 3, bottomCut = 1, bushPreparation = FALSE, filterSOR = FALSE,

                              do.CrownBase = TRUE,
                              groundCutHeight = 0, silent = TRUE, retainPointClouds = FALSE,
                              nr_cores = 0, nr_cores_params = 0,


                              treeNames = NULL,
                              smallTrees = FALSE,
                              dirPath = paste0(getwd(),"/")){



  if(crownParameters){
    writeLAZ <- TRUE # else will get an error
  }

  if(ipad){
    clipHeight <- 2.2
    bottomCut <- 0.2
    selector <- "xyzcitRGB0"
  }
  gstart <- Sys.time()

  fileFinder <- removeUmlaut(fileFinder)

  ### FOLDERS AND FILE NAMES ####
  {
    # DEBUGGING
    if(2==1){
      ipad <- FALSE

      fileFinder = "p14"
      mode = "ALLGO"
      do.CrownBase = TRUE
      voxelSize = 0
      cutWindow = c(-1000,-1000,2000)
      zScale = 1
      limitShare = 0.004

      quantileIntensity = 15
      numberIntensity = 0
      CC_level = 10
      CC_numberPoints = 1000

      clipHeight = 3
      bottomCut = 1
      bushPreparation = FALSE
      filterSOR = FALSE

      silent = TRUE
      dirPath = paste0(getwd(),"/")
      alternativeCrownBase.Ratio <- 0.3



      do.Plot <- T # if F no images are produced

      groundCutHeight <- 0 #how much distance between lowest point of tree and ground (if you cut off the ground)

      # HULL VOLUME PARAMETER
      vol.alpha <- 2 # for rough crowns (balloon)
      #vol.alpha <- 0.4 # OLD SETTING very tight


      # DON'T CHANGE THE OTHER SETTINGS!

      #### SETTINGS FOR TREE HEIGHT MEASUREMENT
      fogFilter.estHeight <- T
      # if TRUE, then height will take 99,99 quantile highest point
      # if FALSE, then height will take maxZ as top z point






    }

    if(!exists("outLAS")){
      outLAS <<- NA # all crowns together in one las file after crownFeeling()
      outLAS_name <<- ""
    }


    if(mode == "COMP"){
      setString <- generateSetString(fileFinder, mode = mode,
                                     threshold = numberIntensity, threshold_percent = quantileIntensity,
                                     bottomCut = bottomCut, clipHeight = clipHeight,
                                     bushPreparation = bushPreparation, filterSOR = filterSOR,
                                     level = CC_level, numberOfPoints = CC_numberPoints, silent = silent)
      dbhPath <- paste0(dirPath,setString,"_dbh/")

      if(!dir.exists(dbhPath)){
        setString <- generateSetString(fileFinder, mode = mode, cutWindow = cutWindow,
                                       threshold = numberIntensity, threshold_percent = quantileIntensity,
                                       bottomCut = bottomCut, clipHeight = clipHeight,
                                       bushPreparation = bushPreparation, filterSOR = filterSOR,
                                       level = CC_level, numberOfPoints = CC_numberPoints, silent = silent)
        dbhPath <- paste0(dirPath,setString,"_dbh/")
      }

    } else
      if(mode == "ALLGO") {
        dbhPath <- generateSetString(fileFinder = fileFinder, mode = mode,
                                     clipHeight = clipHeight, bottomCut = bottomCut,
                                     bushPreparation = bushPreparation,
                                     filterSOR = filterSOR, silent = TRUE)
        dbhPath <- paste0(dirPath,dbhPath,"_dbh/")

        if(!dir.exists(dbhPath)){
          dbhPath <- generateSetString(fileFinder = fileFinder, mode = mode, cutWindow = cutWindow,
                                       clipHeight = clipHeight, bottomCut = bottomCut,
                                       bushPreparation = bushPreparation,
                                       filterSOR = filterSOR, silent = TRUE)
          dbhPath <- paste0(dirPath,dbhPath,"_dbh/")
        }

        if(!dir.exists(dbhPath)){
          warning("There are no clustered input files! \nPlease run functions extractVegetation and clusterSplit/stemSplit before!\n")
          cat("Terminating on",format(Sys.time()),"\n\n")
        }


      }

    crownString <- generateSetString(limitShare = limitShare, zScale = zScale, voxelSize = voxelSize)
    crownPath <- paste0(dbhPath,"crowns",crownString,"/")
    setString_out <- crownString


    locationStr <- generateSetString(cutWindow = cutWindow, silent = silent)
    fileName <- paste0(fileFinder,"_crowns", locationStr,".laz")

    crownPathSingles <- paste0(crownPath,"singleTrees/")
    if(!dir.exists(crownPathSingles)) dir.create(crownPathSingles)


    #crownPathSingles <- paste0(crownPath,"singleTrees",format(Sys.time(), "%Y%m%d_%H%M"),"/") # nice idea timestamp on folder
    crownPathSinglesRef <- paste0(crownPath,"compareQuality/") #not used in this function







     }

  sink(paste0(crownPath,fileFinder,"_computeTreeParams_",format(Sys.time(), "%Y%m%d_%H%M"),"_Rcons.txt"), append = TRUE, split = TRUE)
  cat("COMPUTATION of Single Tree Parameters for case",fileFinder,"\n")
  cat("Job started:",format(gstart, "%Y-%m-%d, %H:%M"),"\n")

  ### READ IN FILES LAS AND TXT ####
  {
    cat("We work with the crown file:",fileName,"\n   in:",crownPath,"\n")
    if(!crownParameters){
      cat("+++FAST ROUND+++ only height estimation and no crown parameters are done.\n")
    }
    cat("\n")

    if(is.na(outLAS) || outLAS_name != paste0(crownPath,fileName)){

      LASfile <- file.path(paste0(crownPath,fileName))

      cat("Reading big crown las file... ")

      rt1 <- Sys.time()
      las <- readLAS(LASfile, select = selector)
      timeNB <- as.difftime(Sys.time() - rt1)
      cat("done in",round(timeNB,1),units(timeNB))
      cat(".\n")


      if(retainPointClouds){
        cat("Retaining outLAS variable for ")
        outLAS <<- las #retain big LAS file in memory, less performant
        outLAS_name <<- paste0(crownPath,fileName)
        cat(outLAS_name,"\n")
      }
    } else {
      cat("We use the prevailing \"outLAS\" variable to extract the crown file...\n")
      las <- outLAS
    }

    trees <- unique(las@data$StemID)
    trees <- sort(trees)

    cat("The", length(trees), "trees range from:",range(trees),"\n")
    beforeLength <- las@header@PHB$`Number of point records`
    cat("Total crown point cloud contains:",thMk(beforeLength),"points.\n")

    # no duplicate filtering at the moment
    #cat("With duplicates there are:",thMk(beforeLength),"points.\n")
    #cat("After filtering...")
    #las <- filter_duplicates(las)
    #afterLength <- las@header@PHB$`Number of point records`
    #cat(" Now left are",thMk(afterLength),"points (removed",thMk(beforeLength - afterLength),"points =",round((beforeLength - afterLength)/beforeLength*100,2),"%).\n")

    stemCoord <- read.table(paste0(dbhPath,"trees_dbh.txt"), header = TRUE)
    stemCoord$dist <- sqrt(stemCoord$x^2 + stemCoord$y^2)


    # cat("We are analyzing the trees from the directory:\n",crownPathSingles,"\n")
    # allTrees <- list.files(crownPathSingles, pattern = ".las", full.names = FALSE)
    # if(length(allTrees) == 0){
    #   cat("No input trees found. \nTerminating.\n\n")
    #   return()
    # }
    cat("\n")
  }

  metaVars <- data.frame("file"=0L, "nPoints"=0,"est.DBH"=0, "est.height"=0, "est.crownLength"=0,
                         "est.crownArea"=0, "est.crownVolume"=0, "est.crownDiameter"=0,
                         "est.DiamCrownMax"=0, "est.DiamCrownMax.Angle"=0,
                         "est.DiamCrownMin"=0, "est.DiamCrownMin.Angle"=0,
                         "est.x.DBH"=0, "est.y.DBH"=0, "est.z.DBH"=0,
                         "est.x.crownDiameter"=0, "est.y.crownDiameter"=0,
                         stringsAsFactors=FALSE) #not implemented yet, was interesting!



  metaVars[1,] = NA
  # metaVarsDBH <- data.frame("file"=0L, "DBH.list"=0, "DBH.rough"=0,
  #                        "DBH.c1"=0, "DBH.c2"=0, "DBH.c3"=0, "DBH.c4"=0, "DBH.c5"=0,
  #                        stringsAsFactors=FALSE) #not implemented yet, was interesting!






  ### INSIDE SETTINGS SECTION ####
  do.Plot <- T # if F no images are produced
  do.Plot <- writePicture


  if(fogFilter.estHeight){
    cat("Fog and snow filter applied to height estimation - working with 99,99 % quantile instead of max Z!\n")
  }




  #cat("These are the trees:\n")
  n.Trees <- length(trees)
  #cat(trees, sep = ", ")
  #cat("\n")

  if(loopEnd == 0 || loopEnd > n.Trees){
    loopEnd <- n.Trees
  }

  if(loopEnd < n.Trees || loopStart != 1){
    cat("\n")
    cat("Subset all trees from", loopStart, "to", loopEnd, "\n")
    cat("The trees to split are:\n   ")
    cat(stemCoord$id[c(loopStart:loopEnd)], sep = ", ")
    cat("\n\n")
  }


  crownBaseList.out <- data.frame("filename" = 0L, "crownBase" = 0)
  crownBase <- NA
  # the criteria for crown base is 1m of more then 2x reference diameter



  gc()
  oSize <- object.size(las)[1] #byte



  if(getRAM){
    totalRAM <- 10000000001 # 10 GB
    tryCatch({
      totalRAM <- benchmarkme::get_ram()[1] #byte
    }, error = function(e) {
      totalRAM <- 10000000001 # 10 GB
    })
    cat("Memory-Info: Size of total pointcloud:", round(oSize/1000^3, 1), "GB, total RAM:",   round(totalRAM/1000^3,1), "GB.\n")
  } else {
    totalRAM <- 10000000001 # 10 GB
  }



  if(nr_cores == 0){
    cat("Automatic number of core setting:\n   ")
    nr_cores <- detectCores()
    cat(nr_cores, "cores available - max ")
    nr_cores <- nr_cores - 1

    maxCores <- floor(totalRAM / oSize / 2) - 1
    cat(maxCores, "cores recommended (RAM ratio of", round(totalRAM / oSize, 1), "/ 2 - 1).\n")
    if(maxCores < 1){
      maxCores <- 1
    }
    if(maxCores <=  nr_cores){
      nr_cores <- maxCores
    }
    cat("\n")
    # if(las@header@PHB$`Number of point records` > 90000000){
    #   cat("LAS File is too big, reducing back to 5 cores!")
    #   if(nr_cores > 5){
    #     nr_cores <- 5
    #   }
    # }
  }

  if(nr_cores > (loopEnd-loopStart+1)){
    nr_cores <- (loopEnd-loopStart+1)
  }

  if(nr_cores == 1){
    # SERIAL SPLITTING

    cat("Running serial now for", (loopEnd-loopStart+1), "trees: \n")
    timePar1 <- Sys.time()
    fdc <<- data.frame()


    for(i in loopStart:loopEnd){


      cat(i,"/",loopEnd,"- tree",trees[i])

                      tstart <- Sys.time()
                      metaVars <- data.frame("file"=0L, "nPoints" = 0, "est.DBH"=0, "est.height"=0,
                                             "est.crownLength"=0, "est.crownBase"=0,
                                             "est.crownArea"=0, "est.crownAreaAlt"=0,
                                             "est.crownVolume"=0, "est.crownDiameter"=0,
                                             "est.DiamCrownMax"=0, "est.DiamCrownMax.Angle"=0,
                                             "est.DiamCrownMin"=0, "est.DiamCrownMin.Angle"=0,
                                             "est.x.DBH"=0, "est.y.DBH"=0, "est.z.DBH"=0,
                                             "est.x.crownDiameter"=0, "est.y.crownDiameter"=0,
                                             "time_secs" = 0,
                                             stringsAsFactors=FALSE) #not implemented yet, was interesting!



                      treeLAS <- filter_poi(las, StemID == trees[i])
                      cat(" split ")
                      tryCatch({

                        {

                          if(is.null(treeNames)){
                            #treeName <- paste0("tree",sprintf("%04d", trees[i]),".laz")
                            treeName <- paste0(sprintf("%04d", trees[i]),"")
                          } else {
                            treeName <- treeNames[i]
                          }


                          #treeName <- paste0("tree",sprintf("%04d", trees[i]),".laz")

                          #### LOADING new tree, SETTINGS ####
                          start <- Sys.time()
                          metaVars$file = treeName


                          nowMeta <- stemCoord[stemCoord$id == trees[i],]

                        }



                        # ESTIMATE HEIGHT
                        #metaVars$est.height <- round(treeLAS@header@PHB$`Max Z` - nowMeta$z + 1.3, 2)

                        up.h <- treeLAS@header@PHB$`Max Z`
                        if(fogFilter.estHeight){
                          up.h <- quantile(treeLAS@data$Z, 0.9999, names = F)
                        }
                        down.h <- nowMeta$z - 1.3
                        diffDraw <- treeLAS@header@PHB$`Max Z` - up.h
                        maxPoint <- treeLAS@data[which(treeLAS@data$Z == treeLAS@header@PHB$`Max Z`),c(1:3)]
                        maxPoint <- maxPoint[1,]
                        metaVars$nPoints <- treeLAS@header@PHB$`Number of point records`
                        metaVars$est.height <- round(up.h - down.h, 2)


                        ## WRITE LAS INDIVIDUAL TREE AND PICTURE ####
                        if(writeLAZ) writeLAS(treeLAS, file = paste0(crownPathSingles, treeName, ".laz"))

                        if(writePicture){
                          if(metaVars$nPoints > 1000){
                            if(smallTrees){
                              treeLAS <- decimate_points(treeLAS, random(10000))
                            } else {
                              treeLAS <- decimate_points(treeLAS, random(4000))
                            }
                          }

                          cat("save ")
                          # Picture of tree from three sides with height measurement
                          png(paste0(crownPathSingles, treeName,".png"), height = 1024, width = 2048, type = "cairo")
                          tryCatch(
                            {
                              par(mfrow=c(1, 3), oma = c(0, 0, 3, 0))
                              plot(treeLAS@data$X, treeLAS@data$Z, asp = 1, cex = 0.001, ylim = c(down.h, up.h))
                              abline(h = c(up.h, down.h), col = "red")
                              lines(x = c(nowMeta$x, maxPoint$X), y = c(down.h, up.h), col = "green", lty = 2)
                              plot(treeLAS@data$Y, treeLAS@data$Z, asp = 1, cex = 0.001, ylim = c(down.h, up.h))
                              abline(h = c(up.h, down.h), col = "red")
                              lines(x = c(nowMeta$y, maxPoint$Y), y = c(down.h, up.h), col = "green", lty = 2)
                              plot(treeLAS@data$X, treeLAS@data$Y, asp = 1, cex = 0.001)
                              lines(x = c(nowMeta$x, maxPoint$X), y = c(nowMeta$y, maxPoint$Y), col = "green", lwd = 2)
                              plotTitle.sin <- paste0(fileFinder, "   ", treeName,
                                                      "   dbh = ", nowMeta$dbh, "cm   est.h = ",
                                                      metaVars$est.height, "m")
                              if(is.element("species", colnames(nowMeta))){
                                if(!is.null(nowMeta$species)){
                                  plotTitle.sin <- paste0(fileFinder, "   ", nowMeta$species,"   ", treeName, "   dbh = ", nowMeta$dbh, "cm   est.h = ",metaVars$est.height, "m")
                                }
                              } else if(is.element("Bart", colnames(nowMeta))){
                                if(!is.null(nowMeta$Bart)){
                                  plotTitle.sin <- paste0(fileFinder, "   ", nowMeta$Bart,"   ", treeName, "   dbh = ", nowMeta$dbh, "cm   est.h = ",metaVars$est.height, "m")
                                }
                              }

                              mtext(plotTitle.sin, outer = TRUE, cex = 2.0)
                              if(drawNumberPoints){
                                mtext(paste0("nPoints = ", thMk(metaVars$nPoints)), outer = TRUE, cex = 1.6, line = -2.7)
                              }
                            }, error=function(e){
                              cat("Error when writing the image, closing device.\n")
                            })
                          dev.off()
                          #par(mfrow=c(1,1))


                        }

                        tstop <- Sys.time()
                        metaVars$time_secs <- round(as.numeric(tstop - tstart, units = "secs"),2)

                        ### LOCAL DENSITY LD100 not implemented anymore ####
                        if(1 == 2){
                          ## Find local density of this tree
                          xCent <- c(nowMeta$x)
                          yCent <- c(nowMeta$y)
                          # trees.info.frame100 <- trees.info[trees.info$X <= xCent + plotSize & trees.info$X >= xCent - plotSize &
                          #                                     trees.info$Y <= yCent + plotSize & trees.info$Y >= yCent - plotSize,]
                          # localDensity100 <- length(trees.info.frame100$ySpan) # trees per hectar in 100 m2 neighborhood
                          # localDensity200 <- length(trees.info.frame$ySpan) # trees per hectar in 100 m2 neighborhood
                          # cat(paste0("\nLocal rectangular densities for ",treeName," is ",localDensity100," trees/100m2 and ",localDensity200," trees/200m2."))

                          # CIRCLE / FIXER PROBEKREIS
                          closest <- nn2(data.frame(xCent, yCent),
                                         query = data.frame(stemCoord$x, stemCoord$y),
                                         k = 1, searchtype = "standard")
                          localDensity100Circle <- length(stemCoord$x[which(closest$nn.dists <= sqrt(100/pi))]) # including the center tree, local density
                          localDensity200Circle <- length(stemCoord$x[which(closest$nn.dists <= sqrt(200/pi))]) # including the center tree, local density
                          #cat(paste0("\nLocal circular densities for ",treeName," is ",localDensity100Circle," trees/100m2 and ",localDensity200Circle," trees/200m2.\n"))
                          metaVars$ld100circ = localDensity100Circle
                          metaVars$ld200circ = localDensity200Circle
                        }




                      }, error = function(e) {
                        # If Error happens at one tree
                        print(paste("MY_ERROR:  ",e))
                        cat("Error catching - closing jpeg device...\n\n")
                        try(dev.off())

                      })


                      cat("done in", metaVars$time_secs, "secs.\n")


                      fdc <<- rbind(fdc, metaVars)
                    }


    metaVars <- fdc
    cat("Serial work done in a ")
    timePar2 <- Sys.time()
    print.difftime(round(timePar2 - timePar1, 1))
    cat("\n")




  } else {
    # PARALLEL SPLITTING



    cl <- makeCluster(nr_cores, outfile="")
    registerDoParallel(cl)
    if(loopEnd == 0){
      loopEnd <- length(trees)
    }


    cat("Going parallel now for", (loopEnd-loopStart+1), "trees on", nr_cores,"cores")
    if(!crownParameters){
      cat(" - ONLY HEIGHTS")
    }
    cat(".\n")
    i <- 1
    timePar1 <- Sys.time()




    #### LOAD BIG LAS AND DIVIDE SINGLE TREES  ####

    fdc <<- foreach(i=loopStart:loopEnd, .errorhandling = 'remove',
                    .packages = c("lidR", "TreeLS", "VoxR"), .verbose = FALSE) %dopar% {



                      tstart <- Sys.time()
                      metaVars <- data.frame("file"=0L, "nPoints" = 0, "est.DBH"=0, "est.height"=0,
                                             "est.crownLength"=0, "est.crownBase"=0,
                                             "est.crownArea"=0, "est.crownAreaAlt"=0,
                                             "est.crownVolume"=0, "est.crownDiameter"=0,
                                             "est.DiamCrownMax"=0, "est.DiamCrownMax.Angle"=0,
                                             "est.DiamCrownMin"=0, "est.DiamCrownMin.Angle"=0,
                                             "est.x.DBH"=0, "est.y.DBH"=0, "est.z.DBH"=0,
                                             "est.x.crownDiameter"=0, "est.y.crownDiameter"=0,
                                             "time_secs" = 0,
                                             stringsAsFactors=FALSE) #not implemented yet, was interesting!

                      #cat("Let's go!\n\n")


                      treeLAS <- filter_poi(las, StemID == trees[i])


                      # # going for every tree taking time
                      # for(i in 1:loopEnd){
                      #   if(i != 1){
                      #     crownBaseList.out <- rbind(crownBaseList.out, data.frame("filename" = allTrees[i-1], "crownBase" = crownBase))
                      #     crownBase <- NA
                      #   }
                      tryCatch({

                        {

                          if(is.null(treeNames)){
                            #treeName <- paste0("tree",sprintf("%04d", trees[i]),".laz")
                            treeName <- paste0(sprintf("%04d", trees[i]),"")
                          } else {
                            treeName <- treeNames[i]
                          }


                          #treeName <- paste0("tree",sprintf("%04d", trees[i]),".laz")

                          #### LOADING new tree, SETTINGS ####
                          start <- Sys.time()
                          metaVars$file = treeName


                          nowMeta <- stemCoord[stemCoord$id == trees[i],]

                        }



                        # ESTIMATE HEIGHT
                        #metaVars$est.height <- round(treeLAS@header@PHB$`Max Z` - nowMeta$z + 1.3, 2)

                        up.h <- treeLAS@header@PHB$`Max Z`
                        if(fogFilter.estHeight){
                          up.h <- quantile(treeLAS@data$Z, 0.9999, names = F)
                        }
                        down.h <- nowMeta$z - 1.3
                        diffDraw <- treeLAS@header@PHB$`Max Z` - up.h
                        maxPoint <- treeLAS@data[which(treeLAS@data$Z == treeLAS@header@PHB$`Max Z`),c(1:3)]
                        maxPoint <- maxPoint[1,]
                        metaVars$nPoints <- treeLAS@header@PHB$`Number of point records`
                        metaVars$est.height <- round(up.h - down.h, 2)


                        ## WRITE LAS INDIVIDUAL TREE AND PICTURE ####
                        if(writeLAZ) writeLAS(treeLAS, file = paste0(crownPathSingles, treeName, ".laz"))

                        if(writePicture){
                          if(metaVars$nPoints > 1000){
                            if(smallTrees){
                              treeLAS <- decimate_points(treeLAS, random(10000))
                            } else {
                              treeLAS <- decimate_points(treeLAS, random(4000))
                            }
                          }
                          # Picture of tree from three sides with height measurement
                          png(paste0(crownPathSingles, treeName,".png"), height = 1024, width = 2048, type = "cairo")
                          tryCatch(
                            {
                              par(mfrow=c(1, 3), oma = c(0, 0, 3, 0))
                              plot(treeLAS@data$X, treeLAS@data$Z, asp = 1, cex = 0.001, ylim = c(down.h, up.h))
                              abline(h = c(up.h, down.h), col = "red")
                              lines(x = c(nowMeta$x, maxPoint$X), y = c(down.h, up.h), col = "green", lty = 2)
                              plot(treeLAS@data$Y, treeLAS@data$Z, asp = 1, cex = 0.001, ylim = c(down.h, up.h))
                              abline(h = c(up.h, down.h), col = "red")
                              lines(x = c(nowMeta$y, maxPoint$Y), y = c(down.h, up.h), col = "green", lty = 2)
                              plot(treeLAS@data$X, treeLAS@data$Y, asp = 1, cex = 0.001)
                              lines(x = c(nowMeta$x, maxPoint$X), y = c(nowMeta$y, maxPoint$Y), col = "green", lwd = 2)
                              plotTitle.sin <- paste0(fileFinder, "   ", treeName, "   dbh = ", round(nowMeta$dbh,2), "cm   est.h = ",metaVars$est.height, "m")
                              if(is.element("species", colnames(nowMeta))){
                                if(!is.null(nowMeta$species)){
                                  plotTitle.sin <- paste0(fileFinder, "   ", nowMeta$species,"   ", treeName, "   dbh = ", round(nowMeta$dbh,2), "cm   est.h = ",metaVars$est.height, "m")
                                }
                              } else if(is.element("Bart", colnames(nowMeta))){
                                if(!is.null(nowMeta$Bart)){
                                  plotTitle.sin <- paste0(fileFinder, "   ", nowMeta$Bart,"   ", treeName, "   dbh = ", round(nowMeta$dbh,2), "cm   est.h = ",metaVars$est.height, "m")
                                }
                              }

                              mtext(plotTitle.sin, outer = TRUE, cex = 2.0)
                              if(drawNumberPoints){
                                mtext(paste0("nPoints = ", thMk(metaVars$nPoints)), outer = TRUE, cex = 1.6, line = -2.7)
                              }

                            }, error=function(e){
                              cat("Error when writing the image, closing device.\n")
                            })
                          dev.off()
                          #par(mfrow=c(1,1))


                        }

                        tstop <- Sys.time()
                        metaVars$time_secs <- round(as.numeric(tstop - tstart, units = "secs"),2)

                        ### LOCAL DENSITY LD100 not implemented anymore ####
                        if(1 == 2){
                          ## Find local density of this tree
                          xCent <- c(nowMeta$x)
                          yCent <- c(nowMeta$y)
                          # trees.info.frame100 <- trees.info[trees.info$X <= xCent + plotSize & trees.info$X >= xCent - plotSize &
                          #                                     trees.info$Y <= yCent + plotSize & trees.info$Y >= yCent - plotSize,]
                          # localDensity100 <- length(trees.info.frame100$ySpan) # trees per hectar in 100 m2 neighborhood
                          # localDensity200 <- length(trees.info.frame$ySpan) # trees per hectar in 100 m2 neighborhood
                          # cat(paste0("\nLocal rectangular densities for ",treeName," is ",localDensity100," trees/100m2 and ",localDensity200," trees/200m2."))

                          # CIRCLE / FIXER PROBEKREIS
                          closest <- nn2(data.frame(xCent, yCent),
                                         query = data.frame(stemCoord$x, stemCoord$y),
                                         k = 1, searchtype = "standard")
                          localDensity100Circle <- length(stemCoord$x[which(closest$nn.dists <= sqrt(100/pi))]) # including the center tree, local density
                          localDensity200Circle <- length(stemCoord$x[which(closest$nn.dists <= sqrt(200/pi))]) # including the center tree, local density
                          #cat(paste0("\nLocal circular densities for ",treeName," is ",localDensity100Circle," trees/100m2 and ",localDensity200Circle," trees/200m2.\n"))
                          metaVars$ld100circ = localDensity100Circle
                          metaVars$ld200circ = localDensity200Circle
                        }




                      }, error = function(e) {
                        # If Error happens at one tree
                        print(paste("MY_ERROR:  ",e))
                        cat("Error catching - closing jpeg device...\n\n")
                        try(dev.off())

                      })
                      return(metaVars)
                    }


    stopCluster(cl)
    getDoParWorkers()

    cat("Parallel work done in a ")
    timePar2 <- Sys.time()
    print.difftime(round(timePar2 - timePar1, 1))
    cat("\n")




    skippedTrees <- sum(lengths(fdc) == 0)
    if(skippedTrees > 0){
      cat("There were", skippedTrees, "trees out of the circle area.\n")

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
      for(j in c(2:length(df[1,]))){ # all num unless tile
        df[,j] <- as.numeric(as.character(df[,j]))
      }
    }

    metaVars <- df


  }





  crownBaseList.out <- data.frame("filename" =  metaVars$file, "crownBase" =  metaVars$est.crownBase)


  if(is.null(treeNames)){
    if(startsWith(metaVars$file[1], "tree")){
      metaVars$id <- as.numeric(strRep(substr(metaVars$file, 5, 99), ".laz", ""))
    } else {
      metaVars$id <- as.numeric(strRep(metaVars$file, ".laz", ""))
    }
  } else {
    metaVars$name <- (metaVars$file)
    metaVars$id <- c(1:length(metaVars$file))
  }
  checkList <- (merge(metaVars, stemCoord, by = "id"))
  write.table(checkList, file = paste0(crownPath, fileFinder, "_metaCrowns_all_",format(Sys.time(), "%Y%m%d_%H%M"),".txt"),
              sep = "\t", na = "", row.names = FALSE)

  # long output format for crown parameters also measured
  # checkList <- checkList[,c(which(colnames(checkList) == "x"),
  #                           which(colnames(checkList) == "y"),
  #                           which(colnames(checkList) == "z"),
  #                           which(colnames(checkList) == "id"),
  #                           which(colnames(checkList) == "name"),
  #                           which(colnames(checkList) == "cluster"),
  #                           which(colnames(checkList) == "randomCol"),
  #                           which(colnames(checkList) == "nPoints"),
  #                           which(colnames(checkList) == "dbh"),
  #                           which(colnames(checkList) == "d.circ"),
  #                           which(colnames(checkList) == "d.circ2"),
  #                           which(colnames(checkList) == "d.ell"),
  #                           which(colnames(checkList) == "d.gam"),
  #                           which(colnames(checkList) == "tegam.d_gam"),
  #                           which(colnames(checkList) == "est.height"),
  #                           which(colnames(checkList) == "est.crownBase"),
  #                           which(colnames(checkList) == "est.crownLength"),
  #                           which(colnames(checkList) == "est.crownArea"),
  #                           which(colnames(checkList) == "est.crownAreaAlt"),
  #                           which(colnames(checkList) == "est.crownVolume"),
  #                           which(colnames(checkList) == "est.crownDiameter"),
  #                           which(colnames(checkList) == "ld100circ"),
  #                           which(colnames(checkList) == "ld200circ"),
  #                           which(colnames(checkList) == "species"),
  #                           which(colnames(checkList) == "inside"),
  #                           which(colnames(checkList) == "time_secs"))
  # ]

    checkList <- checkList[,c(which(colnames(checkList) == "x"),
                              which(colnames(checkList) == "y"),
                              which(colnames(checkList) == "z"),
                              which(colnames(checkList) == "id"),
                              which(colnames(checkList) == "cluster"),
                              which(colnames(checkList) == "randomCol"),
                              which(colnames(checkList) == "nPoints"),
                              which(colnames(checkList) == "dbh"),
                              which(colnames(checkList) == "d.circ"),
                              which(colnames(checkList) == "d.circ2"),
                              which(colnames(checkList) == "d.ell"),
                              which(colnames(checkList) == "d.gam"),
                              which(colnames(checkList) == "tegam.d_gam"),
                              which(colnames(checkList) == "est.height"),
                              which(colnames(checkList) == "species"),
                              which(colnames(checkList) == "inside"),
                              which(colnames(checkList) == "time_secs"))
    ]






  #### OUTPUTTING ####
  cat("Tree splitting done in ")
  timeNB <- as.difftime(Sys.time() - gstart)
  cat(round(timeNB,1),units(timeNB),"\n")


  # if you want to compare different height settings, enable this again
  # write.table(checkList, file = paste0(dbhPath, "trees_height_",setString_out,".txt"),
  #             sep = "\t", na = "", row.names = FALSE)
  write.table(checkList, file = paste0(dbhPath, "trees_height.txt"),
              sep = "\t", na = "", row.names = FALSE)


  #rm(filter.quant, tolerance, start, gstart, loopEnd)
  sink()
  gc()

  if(crownParameters){
    #cat("Zwischenzeit:",format(Sys.time(), "%Y-%m-%d, %H:%M"),"\n\n")

    computeCrownParams(fileFinder = fileFinder, mode = mode, detail.level = detail.level,
                       loopStart = loopStart, loopEnd = loopEnd,
                       maxRadius = maxRadius,ipad = ipad, cutWindow = cutWindow,
                       zScale = zScale, limitShare = limitShare,
                       voxelSize = voxelSize, vol.alpha = vol.alpha,
                       alternativeCrownBase.Ratio = alternativeCrownBase.Ratio,
                       fogFilter.estHeight = fogFilter.estHeight, selector = selector,
                       limitSpanSide = limitSpanSide,
                       limitSpanArea = limitSpanArea, referenceDiameterLimit = referenceDiameterLimit,
                       quantileIntensity = quantileIntensity, numberIntensity = numberIntensity,
                       CC_level = CC_level, CC_numberPoints = CC_numberPoints,
                       clipHeight = clipHeight, bottomCut = bottomCut,
                       bushPreparation = bushPreparation, filterSOR = filterSOR,
                       do.CrownBase = do.CrownBase, groundCutHeight = groundCutHeight,
                       silent = silent, retainPointClouds = retainPointClouds,
                       nr_cores = nr_cores_params, treeNames = treeNames,
                       smallTrees = smallTrees, dirPath = dirPath,
                       drawImages = F)



    gstop <- Sys.time()
    cat("Global ")
    print.difftime(round(gstop - gstart,1))
    globalTimeDiff <<- paste0(round(gstop-gstart, 1), " ", units(gstop - gstart))

    sink()
    gc()
  }

  cat("\n\n\n")




}




#' @export
plotOneTree <- function(treePath = "D:/tree0354.laz", pic = "cb", stemBase = 12.3){
  tempTree <- readLAS(treePath)
  minZ <- tempTree@header@PHB$`Min Z`
  maxZ <- tempTree@header@PHB$`Max Z`
  lZ <- maxZ - minZ
  plot(tempTree@data$Y, tempTree@data$Z, asp = 1, cex = 0.1, main = treePath)
  abline(h = c(minZ,stemBase), col = "red")
  abline(h = c(maxZ,stemBase + 0.07), col = "black")
  mtext(paste0("stem=",round(stemBase-minZ,1),"m"), col = "red", side = 4, adj = 0, cex = 1.3)
  mtext(paste0("crown=",round(maxZ-stemBase,1),"m(",round((maxZ-stemBase)/lZ*100,0),"%)"), side = 4, adj = 1, cex = 1.3)
}
