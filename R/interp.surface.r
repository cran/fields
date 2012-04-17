# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"interp.surface" <- function(obj, loc) {
    
    # obj is a surface or image  object like the list for contour, persp or image.
    # loc a matrix of 2 d locations -- new points to evaluate the surface.
    x <- obj$x
    y <- obj$y
    z <- obj$z
    nx <- length(x)
    ny <- length(y)
    # this clever idea for finding the intermediate coordinates at the new points
    # is from J-O Irisson
    lx <- approx(x, 1:nx, loc[, 1])$y
    ly <- approx(y, 1:ny, loc[, 2])$y
    lx1 <- floor(lx)
    ly1 <- floor(ly)
    # x and y distances between each new point and the closest grid point in the lower left hand corner.
    ex <- lx - lx1
    ey <- ly - ly1
    # fix up weights to handle the case when loc are equal to
    # last grid point.  These have been set to NA above.
    ex[lx1 == nx] <- 1
    ey[ly1 == ny] <- 1
    lx1[lx1 == nx] <- nx - 1
    ly1[ly1 == ny] <- ny - 1
    # bilinear interpolation finds simple weights based on the
    # the four corners of the grid box containing the new
    # points.
    return(z[cbind(lx1, ly1)] * (1 - ex) * (1 - ey) + z[cbind(lx1 + 
        1, ly1)] * ex * (1 - ey) + z[cbind(lx1, ly1 + 1)] * (1 - 
        ex) * ey + z[cbind(lx1 + 1, ly1 + 1)] * ex * ey)
}
