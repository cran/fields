# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"unscale" <- function(x, x.center, x.scale) {
    x <- scale(x, center = FALSE, scale = 1/x.scale)
    x <- scale(x, center = -x.center, scale = FALSE)
    x
}
