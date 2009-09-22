# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"fastTps" <- function(x, Y, m = NULL, p = NULL, theta, 
    lon.lat = FALSE, ...) {
    x <- as.matrix(x)
    d <- ncol(x)
    if (is.null(p)) {
        if (is.null(m)) {
            m <- max(c(2, ceiling(d/2 + 0.1)))
        }
        p <- (2 * m - d)
        if (p <= 0) {
            stop(" m is too small  you must have 2*m -d >0")
        }
    }
    # special arguments to send to the wendland covariance/taper function.
    # see nearest.dist for some explanation of 'method'
    cov.args <- list(k = p, Dist.args = list(method = ifelse(!lon.lat, 
        "euclidean", "greatcircle")))
    Tpscall <- match.call()
    Tpscall$cov.function <- "Compactly supported Wendland function "
    mKrig(x, Y, cov.function = "wendland.cov", m = m, cov.args = cov.args, 
        theta = theta, find.trA = TRUE, ...)
}
