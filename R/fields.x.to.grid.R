# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"fields.x.to.grid" <- function(x, nx = 80, ny = 80, xy = c(1, 2)) {
    if (is.null(x)) {
       stop("Need a an x matrix to determine ranges for grid")
        }
    M <- ncol(x)
    grid.list <- as.list(1:M)
    # add columns names
    names(grid.list) <- dimnames(x)[[2]]
    #     cruise through x dimensions and find medians.
    for (k in 1:M) {
        grid.list[[k]] <- median(x[, k])
    }
    #
    #
    # overwrite with sequences for the two variables of surface
    xr <- range(x[, xy[1]])
    yr <- range(x[, xy[2]])
    grid.list[[xy[1]]] <- seq(xr[1], xr[2], , nx)
    grid.list[[xy[2]]] <- seq(yr[1], yr[2], , ny)
    grid.list
}
