

#' Convert points into polygon for circle section by ordering to angle - not used in package
#' @param pts spatial points object
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

#' Get best cluster for a given number of best clusters to search
#' Credits to the developers of the edci package in R
#' 
#' @param clust cluster object as from oregMclust or circMclust
#' @param nc number of best cluster
#' @param crit which column extracts best clusters, default: "value"
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

#' Delete duplicated clusters via rounding errors or maximization value of zero
#' Credits to the developers of the edci package in R
#' 
#' @param clust numerical matrix where columns are parameters of clusters
#' @param prec decimal places for rounding, positive integer, optional
#' @param ncol number of columns that gives clusters
#' @param dz when TRUE, the clusters with 0 are deleted
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

#' Debugged function to get circular clusters from edci
#' for parameter help please refer to 
#' https://rdrr.io/cran/edci/man/circMclust.html
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


LMCircleFitDebug <- function (XY, ParIni = NA, LambdaIni = 1, epsilon = 1e-06, IterMAX = 50) 
{
  if (length(ParIni) != 3) 
    ParIni <- NA
  if (!is.numeric(ParIni)) 
    ParIni <- NA
  if (sum(is.na(ParIni))==length(ParIni))
    ParIni <- estimateInitialGuessCircle(XY)
  lambda_sqrt = sqrt(LambdaIni)
  Par = ParIni
  Jtmp = CurrentIteration(Par, XY)
  J <- Jtmp$J
  g <- Jtmp$g
  F <- Jtmp$F
  for (iter in 1:IterMAX) {
    while (TRUE) {
      DelPar = mldivide(rbind(J, lambda_sqrt * diag(3)), 
                        rbind(g, 0, 0, 0))
      progress = sqrt(sum(DelPar^2))/(sqrt(sum(Par^2)) + 
                                        epsilon)
      if (progress < epsilon) 
        break
      ParTemp = Par - t(DelPar)
      Jtmp = CurrentIteration(ParTemp, XY)
      JTemp <- Jtmp$J
      gTemp <- Jtmp$g
      FTemp <- Jtmp$F
      if (FTemp < F && ParTemp[3] > 0) {
        lambda_sqrt = lambda_sqrt/2
        break
      }
      else {
        lambda_sqrt = lambda_sqrt * 2
        next
      }
    }
    if (progress < epsilon) 
      break
    Par = ParTemp
    J = JTemp
    g = gTemp
    F = FTemp
  }
  Par
}