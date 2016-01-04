# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"summary.spatial.design" <- function(object, digits = 4, 
    ...) {
    x <- object
    class(x) <- ("summary.spatial.design")
    x
}
