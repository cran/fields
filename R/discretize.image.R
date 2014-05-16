# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"discretize.image" <- function(x, m = 64, n = 64, 
    grid = NULL, expand = c(1, 1), boundary.grid = FALSE, na.rm=TRUE) {
    #
    # set up discretized grid based on x
    #
    out <- list()
       
    if (length(expand) == 1) 
        expand <- rep(expand, 2)
    if (is.null(grid)) {
        grid <- list()
        xr <- range(x[, 1], na.rm = na.rm)
        deltemp <- (xr[2] - xr[1]) * (expand[1] - 1) * 0.5
        grid$x <- seq(xr[1] - deltemp, xr[2] + deltemp, , m)
        yr <- range(x[, 2], na.rm = na.rm)
        deltemp <- (yr[2] - yr[1]) * (expand[2] - 1) * 0.5
        grid$y <- seq(yr[1] - deltemp, yr[2] + deltemp, , n)
    }
    # find cut points for boundaries assuming midpoints
    if (!boundary.grid) {
        xcut <- fields.convert.grid(grid$x)
        ycut <- fields.convert.grid(grid$y)
    }
    else {
        # cut points given boundaries
        xcut <- grid$x
        ycut <- grid$y
    }
    # locate bin ids for each location
    index <- list( as.numeric(cut(x[, 1], xcut)), as.numeric(cut(x[, 2], ycut)))
    m <- length(xcut) - 1
    n <- length(ycut) - 1
    grid <- grid
    # 2 d histogram
    hist<- matrix( 0, m, n)
    tempHist<- table( index[[1]], index[[2]])
    ix<- as.numeric(dimnames( tempHist)[[1]])
    iy<- as.numeric(dimnames( tempHist)[[2]])    
    if (!boundary.grid) {
    #discretized locations
        loc <- cbind( grid$x[ index[[1]] ], grid$y[ index[[2]] ] )  
    }
    else {
        out$loc <- NA
    }
    return( list( m=m,n=n, grid=grid, index=index, ix= ix, iy=iy, hist=hist, loc=loc) )
}
