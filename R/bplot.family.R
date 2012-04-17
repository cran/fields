# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


"bplot.xy" <- function(x, y, N = 10, breaks = pretty(x, 
    N, eps.correct = 1), plot = TRUE, ...) {
    NBIN <- length(breaks) - 1
    centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
    obj <- split(y, cut(x, breaks))
    if (length(obj) == 0) {
        stop("No points within breaks")
    }
    if (plot) {
        bplot(obj, at = centers, show.names = FALSE, axes = TRUE, 
            ...)
        axis(1)
    }
    else {
        return(list(centers = centers, breaks = breaks, boxplot.obj = boxplot(obj, 
            plot = FALSE)))
    }
}

bplot <- function(x, by, pos = NULL, at = pos, add = FALSE, 
    boxwex = 0.8, xlim = NULL, ...) {
    if (!missing(by)) {
        x <- split(c(x), as.factor(by))
    }
    if (!add & !is.null(at) & is.null(xlim)) {
        xlim <- range(at)
    }
    if (!is.null(at)) {
        boxwex <- boxwex * min(diff(sort(at)))
    }
    boxplot(x, at = at, xlim = xlim, add = add, boxwex = boxwex, 
        ...)
    
}


