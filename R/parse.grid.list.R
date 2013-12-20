# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"parse.grid.list" <- function(grid.list, order.variables = "xy") {
    #
    # utility to find the x and y sequences in grid.list
    # this is used in predictSurface and as.surface
    #
    M <- length(grid.list)
    gcounts <- unlist(lapply(grid.list, FUN = length))
    xy <- (1:M)[gcounts > 1]
    if (length(xy) > 2) {
        stop("only two components of the grid list\ncan have more than one element")
    }
    #
    # swap the roles of x and y
    if (order.variables == "yx") {
        xy <- xy[2:1]
    }
    #
    #
    # here is the good stuff
    #
    nx <- gcounts[xy[1]]
    ny <- gcounts[xy[2]]
    x <- grid.list[[xy[1]]]
    y <- grid.list[[xy[2]]]
    #
    #  extract the names of the x and y components of the
    #  list
    #
    xlab <- names(grid.list)[xy[1]]
    ylab <- names(grid.list)[xy[2]]
    xlab <- ifelse(is.null(xlab), "X", xlab)
    ylab <- ifelse(is.null(ylab), "Y", ylab)
    list(x = x, y = y, nx = nx, ny = ny, xlab = xlab, ylab = ylab, 
        xy = xy)
}
