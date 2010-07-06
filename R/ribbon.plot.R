# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
ribbon.plot <- function(x, y, z, zlim = NULL, col = tim.colors(256), 
    transparent.color = "white", ...) {
    N <- length(x)
    x1 <- (x[1:(N - 1)] + x[2:(N)])/2
    y1 <- (y[1:(N - 1)] + y[2:(N)])/2
    x1 <- c(x[1] - (x[2] - x[1])/2, x1, x[N] + (x[N] - x[N - 
        1])/2)
    y1 <- c(y[1] - (y[2] - y[1])/2, y1, y[N] + (y[N] - y[N - 
        1])/2)
    eps <- 1e-07
    if (is.null(zlim)) {
        zlim <- range(c(z), na.rm = TRUE)
    }

# convert z values to a color scale.
    colz <- color.scale(z, col=col,transparent.color=transparent.color)

    segments(x1[1:(N)], y1[1:(N)], x1[2:(N + 1)], y1[2:(N + 1)], 
        col = colz, ...)
}
