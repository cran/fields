# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


"bplot.xy" <- function(x, y, N = 10, breaks = pretty(x, N, eps.correct = 1),
                       plot = TRUE, ...) {
    NBIN <- length(breaks) - 1
    centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
    obj<- split( y, cut(x,breaks))
    if( length( obj)==0){ stop("No points within breaks")}
    if (plot) {
      bplot(obj, at = centers,show.names=FALSE,axes=TRUE,...)
      axis(1)
    }
    else{
      return(list(centers=centers, breaks=breaks,
                   boxplot.obj=boxplot(obj, plot=FALSE)))}
}

bplot<- function(x, by,pos=NULL, at= pos, add=FALSE,boxwex=.8, ...){
   if (!missing(by)) { 
     x <- split(c(x), as.factor(by))}
   if( !add & !is.null(at)){
     xlim<- range(at)}
   else{
     xlim<- NULL}
   if( !is.null(at)){
     boxwex<- boxwex* min(diff( sort( at))) }
   
   boxplot( x,at=at,xlim=xlim, add=add,boxwex=boxwex,...)
   
  }


########################### these are old bplot utilities that
########################### may be useful but are not needed for
########################### current bplot and bplot.xy
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#"bplot" <- function(x, by, style = "tukey", outlier = TRUE, 
#    plot = TRUE, ...) {
#    obj <- stats.bplot(x, style = style, outlier = outlier, by = by)
#    if (plot) {
#        bplot.obj(obj, outlier = outlier, ...)
#    }
#    else {
#        return(obj)
#    }
#    invisible()
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
#"bplot.xy" <- function(x, y, N = 10, breaks = pretty(x, 
#    N, eps.correct = 1), style = "tukey", outlier = TRUE, plot = TRUE, 
#    xaxt = "s", ...) {
#    out <- list()
#    NBIN <- length(breaks) - 1
#    centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
#    obj <- as.list(1:NBIN)
#    names(obj) <- format(1:NBIN)
#    for (k in 1:NBIN) {
#        obj[[k]] <- describe.bplot(y[x < breaks[k + 1] & x > 
#            breaks[k]], style = style, outlier = outlier)
#    }
#    if (plot) {
#        bplot.obj(obj, pos = centers, label.cex = 0, outlier = outlier, 
#            , xaxt = xaxt, ...)
#    }
#    else {
#        return(list(centers = centers, breaks = breaks, bplot.obj = obj))
#    }
#    invisible()
#}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"draw.bplot" <- function(temp, width, xpos, outlier = TRUE, 
    style = "tukey") {
    if (temp$N < 1) 
        return()
    if (style == "quantile") {
        temp <- temp[!is.na(temp)]
        quant <- c(0.05, 0.25, 0.5, 0.75, 0.95)
        bb <- quantile(temp, quant)
        mid <- xpos
        low <- mid - width * 0.5
        high <- mid + width * 0.5
        if (length(temp) > 5) {
            y <- c(bb[1], bb[1], NA, bb[1], bb[2], NA, bb[2], 
                bb[2], bb[4])
            x <- c(high, low, NA, mid, mid, NA, high, low, low)
            y <- c(y, bb[4], bb[2], bb[3], bb[3], NA, bb[4], 
                bb[5], bb[5], bb[5])
            x <- c(x, high, high, high, low, NA, mid, mid, high, 
                low)
            lines(x, y)
        }
        if (length(temp) > 5) {
            outs <- temp[(temp < bb[1]) | (temp > bb[5])]
        }
        else outs <- temp
        olen <- length(outs)
        if ((olen > 0) & outlier) 
            points(rep(mid, olen), outs)
    }
    if (style == "tukey") {
        temp <- temp[!is.na(temp)]
        quant <- c(0.05, 0.25, 0.5, 0.75, 0.95)
        bb <- quantile(temp, quant)
        iqr <- bb[4] - bb[2]
        mid <- xpos
        low <- mid - width * 0.5
        high <- mid + width * 0.5
        bb[1] <- min(temp[temp >= bb[2] - 1.5 * iqr])
        bb[5] <- max(temp[temp <= bb[4] + 1.5 * iqr])
        if (length(temp) > 5) {
            y <- c(bb[1], bb[1], NA, bb[1], bb[2], NA, bb[2], 
                bb[2], bb[4])
            x <- c(high, low, NA, mid, mid, NA, high, low, low)
            y <- c(y, bb[4], bb[2], bb[3], bb[3], NA, bb[4], 
                bb[5], bb[5], bb[5])
            x <- c(x, high, high, high, low, NA, mid, mid, high, 
                low)
            lines(x, y)
        }
        if (length(temp) > 5) {
            outs <- temp[(temp < bb[2] - 3 * iqr) | (temp > bb[4] + 
                3 * iqr)]
        }
        else outs <- temp
        olen <- length(outs)
        if ((olen > 0) & outlier) 
            points(rep(mid, olen), outs)
    }
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"draw.bplot.obj" <- function(obj, width, xpos, outlier = TRUE, 
    horizontal = FALSE, lwd = NA, col = NA) {
    #
    # fill in defaults if not specfied
    if (missing(lwd)) 
        lwd <- par()$lwd
    if (missing(col)) 
        lwd <- par()$col
    N <- obj$N
    bb <- obj$bb
    mid <- xpos
    low <- mid - width * 0.5
    high <- mid + width * 0.5
    if (N > 5) {
        y <- c(bb[1], bb[1], NA, bb[1], bb[2], NA, bb[2], bb[2], 
            bb[4])
        x <- c(high, low, NA, mid, mid, NA, high, low, low)
        y <- c(y, bb[4], bb[2], bb[3], bb[3], NA, bb[4], bb[5], 
            bb[5], bb[5])
        x <- c(x, high, high, high, low, NA, mid, mid, high, 
            low)
        if (horizontal) {
            lines(y, x, lwd = lwd, col = col)
        }
        else {
            lines(x, y, lwd = lwd, col = col)
        }
    }
    outs <- obj$out
    olen <- length(outs)
    if ((olen > 0) & outlier) {
        if (horizontal) {
            points(outs, rep(mid, olen), col = col)
        }
        else {
            points(rep(mid, olen), outs, col = col)
        }
    }
}



# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"bplot.obj" <- function(data, pos = NA, width = NULL, 
    labels = NULL, las = NULL, add = FALSE, space = 0.25, sort.names = FALSE, 
    xlab = "", ylab = "", label.cex = 1, xaxt = "n", outlier = TRUE, 
    horizontal = FALSE, lwd = NA, col = NA, ...) {
    cols <- length(data)
    # figure out what colors and line widths to use for bplots
    # if fewer of these than cols are given then rep them to fill
    # out a vector of length cols
    if ((is.na(col))[1]) {
        col <- par()$col
    }
    if (length(col) < cols) {
        col <- (rep(col, cols))[1:cols]
    }
    if (is.na(lwd)) {
        lwd <- par()$lwd
    }
    if (length(lwd) < cols) {
        lwd <- (rep(lwd, cols))[1:cols]
    }
    range.data <- c(NA, NA)
    if (is.null(labels)) {
        labels <- names(data)
    }
    if (is.na(pos[1])) {
        pos <- 1:cols
        if (sort.names) {
            pos <- order(labels)
        }
    }
    if (is.null(width)) {
        width <- min(diff(sort(pos))) * space
        if (cols == 1) 
            width <- space
    }
    if (length(width) == 1) 
        width <- rep(width, cols)
    # determine limits for plotting bplots
    # and set up plot if bplots are not being added to an existing plot
    if (!add) {
        for (k in 1:cols) {
            # if plotting outliers too use the full range
            # it outliers are not plotted then just use limits of boxplot whiskers (in $bb)
            if (outlier) {
                range.data <- range(c(range.data, data[[k]]$range), 
                  na.rm = TRUE)
            }
            else {
                range.data <- range(c(range.data, data[[k]]$bb), 
                  na.rm = TRUE)
            }
        }
        temp1 <- range.data
        temp2 <- range(c(pos - (0.5 * width)/space, pos + (0.5 * 
            width)/space))
        # set up empty plot
        if (horizontal) {
            plot(temp1, temp2, type = "n", yaxt = xaxt, xlab = xlab, 
                ylab = ylab, ...)
        }
        else {
            plot(temp2, temp1, type = "n", xaxt = xaxt, xlab = xlab, 
                ylab = ylab, ...)
        }
    }
    # loop through data and draw each bplot
    for (i in 1:cols) {
        draw.bplot.obj(data[[i]], width[i], pos[i], outlier = outlier, 
            horizontal = horizontal, lwd = lwd[i], col = col[i])
    }
    # add labels
    if (label.cex > 0) {
        if (is.null(las)) {
            if (length(labels) > 7) {
                las <- 2
            }
            else {
                las <- 1
            }
        }
        if (horizontal) {
            axis.loc <- 2
        }
        else {
            axis.loc <- 1
        }
        axis(axis.loc, pos, labels, tick = FALSE, las = las, 
            adj = 0.5, cex = label.cex)
    }
    invisible()
}

