# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
poly.image.regrid <- function(x) {
    ##################################
    temp.addcol <- function(X) {
        N <- ncol(X)
        # add extra columns to X, on either side
        cbind(X[, 1] - (X[, 2] - X[, 1]), X, (X[, N] - X[, (N - 
            1)]) + X[, N])
    }
    ###############################
    #     find approximate grid with z values at centers
    M <- nrow(x)
    N <- ncol(x)
    # new x matrix that is the midpoints of original grid points.
    x <- (x[, 1:(N - 1)] + x[, 2:N])/2
    x <- (x[1:(M - 1), ] + x[2:M, ])/2
    # now add extra rows and columns on all sides
    x <- t(temp.addcol(x))
    t(temp.addcol(x))
}
poly.image <- function(x, y, z, col = tim.colors(64), 
    transparent.color = "white", midpoint = FALSE, zlim = range(z, 
        na.rm = TRUE), xlim = range(x), ylim = range(y), add = FALSE, 
    border = NA, ...) {
    # check dimensions
    Dx <- dim(x)
    Dy <- dim(y)
    if (any((Dx - Dy) != 0)) {
        stop(" x and y matrices should have same dimensions")
    }
    # check whether grid and z values coincide.
    Dz <- dim(z)
    if (all((Dx - Dz) == 0) & !midpoint) {
        # expand grid in a linear way so that the z are not
        # grid box centers
        x <- poly.image.regrid(x)
        y <- poly.image.regrid(y)
    }
    # code values in z based on range to colors.
    # if midpoint is true z values will be averaged first
    zcol <- drape.color(z, col = col, midpoint = midpoint, zlim = zlim, 
        transparent.color = transparent.color)
    # blank if not adding to an exising plot
    if (!add) {
        plot(xlim, ylim, type = "n", ...)
    }
    N <- ncol(x)
    Nm1 <- N - 1
    M <- nrow(x)
    Mm1 <- M - 1
    for (i in (1:Mm1)) {
        # draw polygons one row at a time
        # uses feature of polygon to draw multiple regions with NAs
        # marking start and end.
        xp <- cbind(x[i, 1:Nm1], x[i + 1, 1:Nm1], x[i + 1, 2:N], 
            x[i, 2:N], rep(NA, Nm1))
        yp <- cbind(y[i, 1:Nm1], y[i + 1, 1:Nm1], y[i + 1, 2:N], 
            y[i, 2:N], rep(NA, Nm1))
        xp <- c(t(xp))
        yp <- c(t(yp))
        polygon(xp, yp, border = NA, col = c(zcol[i, 1:Nm1]))
    }
}