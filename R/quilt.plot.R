# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"quilt.plot" <- function(x, y, z, nx=64, ny=64,nrow = nx, ncol = ny, 
    grid = NULL, add.legend = TRUE, add = FALSE, col = tim.colors(256), 
    ...) {
#
# note that nrow and ncol refer to the resulting "image format" for plotting. 
# here the x values are the rows and the y values are the columns
# 
    x <- as.matrix(x)
    if (ncol(x) == 2) {
        z <- y
    }
    if (ncol(x) == 1) {
        x <- cbind(x, y)
    }
    if (ncol(x) == 3) {
        z <- x[, 3]
        x <- x[, 1:2]
    }
    # at this point x should be a 2 column matrix of x-y locations
    #  z is a vector or one column matrix of the z values.
    #discretize data
    out.p <- as.image(z, x = x, nrow = nrow, ncol = ncol,
                            na.rm = TRUE, grid=grid)
    #plot it
    if (add.legend) {
        image.plot(out.p, col = col, add = add, ...)
    }
    else {
        image(out.p, col = col, add = add, ...)
    }
}
