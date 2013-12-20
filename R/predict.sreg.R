# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"predict.sreg" <- function(object, x, derivative = 0, 
    model = 1, ...) {
    if (missing(x)) {
        x <- object$x
    }
    c(splint(object$predicted$x, object$predicted$y[, model], 
        x, derivative = derivative, ...))
}
