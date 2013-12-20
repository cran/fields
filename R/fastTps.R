# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"fastTps" <- function(x, Y, m = NULL, p = NULL, theta, 
    lon.lat = FALSE, find.trA=TRUE,lambda=0,...) {
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
    if( lambda==0){
      warning("fastTps will interpolate observations")}
    mKrig(x, Y, cov.function = "wendland.cov", m = m, cov.args = cov.args, 
        theta = theta, find.trA = find.trA,lambda=lambda, ...)
}
