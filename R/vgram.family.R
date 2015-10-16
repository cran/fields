"vgram" <- function(loc, y, id = NULL, d = NULL, lon.lat = FALSE, 
                    dmax = NULL, N = NULL, breaks = NULL, 
                    type=c("variogram", "covariogram", "correlogram")) {
  
  type=match.arg(type)
  
  # coerce to matrix
  y <- cbind(y)
  # if nearest neighbor indices are missing create all possible pairs.
  if (is.null(id)) {
    n <- nrow(loc)
    is = rep(1:n, n)
    js = rep(1:n, rep(n, n))
    ind <- is > js
    id <- cbind(is, js)[ind, ]
  }
  
  # if distances are missing calculate these
  if (is.null(d)) {
    loc <- as.matrix(loc)
    if (lon.lat) {
      d <- rdist.earth.vec(loc[id[,1],], loc[id[,2],]) #we want result in miles, not meters
    }
    else {
      d <- rdist.vec(loc[id[,1],], loc[id[,2],])
    }
  }
  
  # normalize columns to create correlogram, if necessary
  #
  if(type == "correlogram") {
    sigma = apply(y, 2, sd, na.rm=TRUE)
    y = sweep(y, 2, (1/sigma), FUN="*")
  }
  
  # center the columns by their mean and get row means if y is a matrix
  #
  colMeans <- apply(y, 2, mean, na.rm=TRUE)
  yCntr = sweep(y, 2, colMeans) 
  y1Cntr = yCntr[id[,1],]
  y2Cntr = yCntr[id[,2],]
  if(type == "variogram") {
    vg <- 0.5 * rowMeans(cbind((y1Cntr - y2Cntr)^2), 
                         na.rm = TRUE)
  }
  else {
    vg <- rowMeans(cbind(y1Cntr * y2Cntr), 
                   na.rm = TRUE)
  }
  #
  #information for returned object
  #
  call <- match.call()
  if (is.null(dmax)) {
    dmax <- max(d)
  }
  od <- order(d)
  d <- d[od]
  vg <- vg[od]
  ind <- d <= dmax & !is.na(vg)
  
  ## add a binned  variogram if breaks are supplied
  out <- list(d = d[ind], vgram = vg[ind], call = call, type=type)
  if (!is.null(breaks) | !is.null(N)) {
    out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
  }
  class(out) = c("vgram", class(out))
  out
}

#calculating cross-covariogram and cross-correlogram (cross-covariance and 
#cross-correlation)
crossCoVGram = function(loc1, loc2, y1, y2, id = NULL, d = NULL, lon.lat = FALSE, 
                        dmax = NULL, N = NULL, breaks = NULL, 
                        type=c("cross-covariogram", "cross-correlogram")) {
  
  type=match.arg(type)
  
  # coerce to matrix
  y1 <- cbind(y1)
  y2 <- cbind(y2)
  
  # if nearest neighbor indices are missing create all possible pairs.
  if (is.null(id)) {
    n1 <- nrow(loc1)
    n2 <- nrow(loc2)
    id <- cbind(rep(1:n1, n2), rep(1:n2, rep(n1, n2)))
  }
  
  # if distances are missing calculate these
  if (is.null(d)) {
    loc1 <- as.matrix(loc1)
    loc2 <- as.matrix(loc2)
    if (lon.lat) {
      d <- rdist.earth.vec(loc1[id[,1],], loc2[id[,2],]) #we want result in miles, not meters
    }
    else {
      d <- rdist.vec(loc1[id[,1],], loc2[id[,2],])
    }
  }
  #
  # calculating covariogram will center the columns by their mean and get row means if y is a matrix
  #
  colMeans1 <- apply(y1, 2, mean, na.rm=TRUE)
  colMeans2 <- apply(y2, 2, mean, na.rm=TRUE)
  y1Cntr = sweep(data.matrix(y1), 2, colMeans1)  # subtract the column means
  y2Cntr = sweep(data.matrix(y2), 2, colMeans2)  # subtract the column means
  #
  # normalize to create cross-correlogram, if necessary
  #
  if(type == "cross-correlogram") {
    sigma1 = apply(y1Cntr, 2, sd, na.rm=TRUE)
    sigma2 = apply(y2Cntr, 2, sd, na.rm=TRUE)
    y1Cntr = sweep(y1Cntr, 2, 1/sigma1, FUN="*")
    y2Cntr = sweep(y2Cntr, 2, 1/sigma2, FUN="*")
  }
  #
  # calculate covariance for the given points
  #
  y1Cntr = y1Cntr[id[,1],]
  y2Cntr = y2Cntr[id[,2],]
  vg <- rowMeans(cbind(y1Cntr*y2Cntr), na.rm = TRUE)
  #
  #information for returned object
  #
  call <- match.call()
  if (is.null(dmax)) {
    dmax <- max(d)
  }
  od <- order(d)
  d <- d[od]
  vg <- vg[od]
  ind <- d <= dmax & !is.na(vg)
  ## add a binned  variogram if breaks are supplied
  out <- list(d = d[ind], vgram = vg[ind], call = call, type=type)
  if (!is.null(breaks) | !is.null(N)) {
    out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
  }
  class(out) = c("vgram", class(out))
  out
}

#plot only the line of the empirical variogram, where the y coordinates of the line are 
#at the means of the bins
plot.vgram = function(x, N=10, breaks = pretty(x$d, N, eps.correct = 1), add=FALSE, ...) {
  otherArgs = list(...)
  type=x$type
  
  #set y axis label if not set by user
  if(is.null(otherArgs$ylab)) {
    if(type=="variogram")
      ylab = "sqrt(Variance)"
    else if(type == "covariogram" || type=="cross-covariogram")
      ylab = "Covariance"
    else if(type == "correlogram" || type=="cross-correlogram")
      ylab = "Correlation"
    else
      stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  } else {
    ylab = otherArgs$ylab
    otherArgs$ylab = NULL
  }
  
  #set x axis label if not set by user
  if(is.null(otherArgs$xlab))
    xlab = "Distance (Miles)"
  else {
    xlab = otherArgs$xlab
    otherArgs$xlab = NULL
  }
  
  #set plot title if not set by user
  if(is.null(otherArgs$main)) {
    if(type=="variogram")
      main = "Empirical Variogram"
    else if(type=="covariogram")
      main = "Empirical Covariogram"
    else if(type=="correlogram")
      main = "Empirical Correlogram"
    else if(type=="cross-covariogram")
      main = "Empirical Cross-Covariogram"
    else if(type=="cross-correlogram")
      main = "Empirical Cross-Correlogram"
    else
      stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  }
  else {
    main = otherArgs$main
    otherArgs$main = NULL
  }
  
  #set ylim for correlogram if not set by user
  if(is.null(otherArgs$ylim)) {
    if(type == "correlogram" || type=="cross-correlogram")
      ylim = c(-1, 1)
    else
      ylim = NULL
  }
  else {
    ylim = otherArgs$ylim
    otherArgs$ylim = NULL
  }
  
  #set line type if not set by user
  if(is.null(otherArgs$type))
    type = "o"
  else {
    type = otherArgs$type
    otherArgs$type = NULL
  }
  
  #get boxplot object
  dat = do.call("boxplotVGram", c(list(x=x, N=N, breaks=breaks, plot=FALSE), otherArgs))
  
  #get bin centers versus bin means
  centers = dat$centers
  ys = dat$boxplot.obj$stats[3,]
  
  #remove NAs
  notNas = !is.na(ys)
  centers = centers[notNas]
  ys = ys[notNas]
  
  #plot
  if(!add)
    do.call(plot, c(list(centers, ys, main=main, xlab=xlab, ylab=ylab, type=type, ylim=ylim), otherArgs))
  else
    do.call(lines, c(list(centers, ys, main=main, xlab=xlab, ylab=ylab, type=type, ylim=ylim), otherArgs))
}

"boxplotVGram" = function(x, N=10, breaks = pretty(x$d, N, eps.correct = 1), plot=TRUE, ...) {
  dists = x$d
  type=x$type
  if(type == "variogram")
    y = sqrt(x$vgram)
  else
    y = x$vgram
  otherArgs = list(...)
  
  #set y axis label if not set by user
  if(is.null(otherArgs$ylab)) {
    if(type=="variogram")
      ylab = "sqrt(Variance)"
    else if(type == "covariogram" || type=="cross-covariogram")
      ylab = "Covariance"
    else if(type == "correlogram" || type=="cross-correlogram")
      ylab = "Correlation"
    else
      stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  } else {
    ylab = otherArgs$ylab
    otherArgs$ylab = NULL
  }
  
  #set x axis label if not set by user
  if(is.null(otherArgs$xlab))
    xlab = "Distance (Miles)"
  else {
    xlab = otherArgs$xlab
    otherArgs$xlab = NULL
  }
  
  #set plot title if not set by user
  if(is.null(otherArgs$main)) {
    if(type=="variogram")
      main = "Empirical Variogram"
    else if(type=="covariogram")
      main = "Empirical Covariogram"
    else if(type=="correlogram")
      main = "Empirical Correlogram"
    else if(type=="cross-covariogram")
      main = "Empirical Cross-Covariogram"
    else if(type=="cross-correlogram")
      main = "Empirical Cross-Correlogram"
    else
      stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  }
  else {
    main = otherArgs$main
    otherArgs$main = NULL
  }
  
  #set ylim for correlogram if not set by user
  if(is.null(otherArgs$ylim)) {
    if(type == "correlogram" || type=="cross-correlogram")
      ylim = c(-1, 1)
    else
      ylim = NULL
  }
  else {
    ylim = otherArgs$ylim
    otherArgs$ylim = NULL
  }
  
  do.call("bplot.xy", c(list(x=dists, y=y, N=N, breaks=breaks, plot=plot, ylab=ylab, 
                             xlab=xlab, main=main, ylim=ylim), otherArgs))
}
