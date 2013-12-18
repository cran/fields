# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"predict.surface" <- function(object, ...) {
    UseMethod("predict.surface")
}

predict.surface.default<- function(object,...){
   cat("predict.surface is now the function predictSurface")
 }

"predictSurface" <- function(object, grid.list = NA, 
    extrap = FALSE, chull.mask = NA, nx = 80, ny = 80, xy = c(1, 
        2), order.variables = "xy", verbose = FALSE, ...) {
    # without grid.list
    # default is 80X80 grid on frist two variables
    # rest are set to median value of x.
    if (is.na(grid.list)[1]) {
        if (is.null(object$x)) {
            stop("Need a an X matrix in the output object")
        }
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = xy)
    }
    #
    # simple check that grid list matches data
    if (!is.null(object$x)) {
        if (ncol(object$x) != length(grid.list)) {
            stop("grid.list does\nnot match x")
        }
    }
    #
    #
    if (verbose) {
        cat("grid.list")
        print(grid.list)
    }
    # find out where the intended x y sequences are in grid.list
    # and their  respective lengths.
    hold <- parse.grid.list(grid.list, order.variables = order.variables)
    xy <- hold$xy
    nx <- hold$nx
    ny <- hold$ny
    if (verbose) {
        print(xy)
        print(c(nx, ny))
    }
    # output  'z' matrix
    out <- matrix(NA, nx, ny)
    #
    # explicitly loop through  grid row by row to reduce memory
    #
    #    note: a direct method would be:
    #            xg<- make.surface.grid(grid.list)
    #            out <-  as.surface( grid.list,predict.(object, xg, ...))
    #
    gtemp <- grid.list
    # this the position of the x grid in the list
    xloc <- xy[1]
    # fill out row by row
    for (i in 1:nx) {
        #      cat(i, ' ')
        # temporary grid.list with one value in 'x' and all the 'y' grid values
        gtemp[[xloc]] <- grid.list[[xloc]][i]
        # here is the heavy lifting
        out[i, ] <- predict(object, make.surface.grid(gtemp), 
            ...)
    }
    #
    # overwrite out object in usual list format for contour, persp, image plotting
    out <- list(x = grid.list[[xy[1]]], y = grid.list[[xy[2]]], 
        z = out)
    #
    # if extrapolate is FALSE set all values outside convex hull to NA
    if (!extrap) {
        #
        #  all locations on grid
        xg <- make.surface.grid(grid.list)
        if (is.na(chull.mask)) {
            chull.mask <- unique.matrix(object$x[, xy])
        }
        out$z[!in.poly(xg[, xy], xp = chull.mask, convex.hull = TRUE)] <- NA
    }
    #
    out
}
