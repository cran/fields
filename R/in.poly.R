# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"in.poly" <- function(xd, xp, convex.hull = FALSE, 
    inflation = 1e-07) {
    if (convex.hull) {
        xp <- xp[chull(xp), ]
    }
    nd <- as.integer(nrow(xd))
    np <- as.integer(nrow(xp))
    #
    # inflate convex hull slightly to include any points actually on the hull
    #
    if (convex.hull) {
        xm <- matrix(c(mean(xp[, 1]), mean(xp[, 2])), nrow = np, 
            ncol = 2, byrow = TRUE)
        xp <- (xp - xm) * (1 + inflation) + xm
    }
    # Note: inpoly FORTRAN has built in quick reject check to be inside
    # the bounding rectangle of the polygon.
    ind <- .Fortran("inpoly", nd = as.integer(nd), as.single(xd[, 
        1]), as.single(xd[, 2]), np = np, as.single(xp[, 1]), 
        as.single(xp[, 2]), ind = as.integer(rep(-1, nd)))$ind
    as.logical(ind)
}
in.poly.grid <- function(grid.list, xp, convex.hull = FALSE, 
    inflation = 1e-07) {
    # loop through rows of grid to fill out a logical matrix of
    # being in (TRUE) or out (FALSE)
    #
    # this is just to avoid creating the full set of image locations.
    if (convex.hull) {
        xp <- xp[chull(xp), ]
    }
    nx <- length(grid.list$x)
    ny <- length(grid.list$y)
    np <- as.integer(nrow(xp))
    #
    # inflate convex hull slightly to include any points actually on the hull
    #
    if (convex.hull) {
        xm <- matrix(c(mean(xp[, 1]), mean(xp[, 2])), nrow = np, 
            ncol = 2, byrow = TRUE)
        xp <- (xp - xm) * (1 + inflation) + xm
    }
    # Note: inpoly FORTRAN has built in quick reject check to be inside
    # the bounding rectangle of the polygon.
    ind <- .Fortran("igpoly", nx = as.integer(nx), xg = as.single(grid.list$x), 
        ny = as.integer(ny), yg = as.single(grid.list$y), np = np, 
        as.single(xp[, 1]), as.single(xp[, 2]), ind = as.integer(rep(-1, 
            nx * ny)))$ind
    return(matrix(as.logical(ind), nrow = nx, ncol = ny))
}
