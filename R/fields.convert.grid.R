# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"fields.convert.grid" <- function(midpoint.grid) {
    # converts from midpoints of a grid to boundaries
    # x are midpoints of grid
    # this will handle unequally spaced points
    x <- sort(midpoint.grid)
    n <- length(x)
    # interior boundaries
    xi <- (x[2:n] + x[1:(n - 1)])/2
    # first and last.
    x1 <- x[1] - (x[2] - x[1])/2
    xnp1 <- x[n] + (x[n] - x[(n - 1)])/2
    #here you have it ...
    c(x1, xi, xnp1)
}
