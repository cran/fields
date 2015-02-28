# fields, Tools for spatial data
# Copyright 2004-2014, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"Tps" <- function(x, Y, m = NULL, p = NULL, scale.type = "range", 
    lon.lat = FALSE, miles = TRUE, method="GCV", GCV=TRUE, ...) {
    x <- as.matrix(x)
    d <- ncol(x)
    if (is.null(p)) {
        if (is.null(m)) {
            m <- max(c(2, ceiling(d/2 + 0.1)))
        }
        p <- (2 * m - d)
        if (p <= 0) {
            stop(" m is too small  you must have 2*m - dimension >0")
        }
    }
#    Tpscall <- match.call()
    if (!lon.lat) {
#        Tpscall$cov.function <- "Thin plate spline radial basis functions (Rad.cov) "
        obj<- Krig(x, Y, cov.function = Rad.cov, m = m, scale.type = scale.type, 
            p = p, method=method, GCV = GCV, ...)
    }
    else {
        # a different coding of the radial basis functions to use great circle distance.
#        Tpscall$cov.function <- "Thin plate spline radial basis functions (RadialBasis.cov) using great circle distance "
       obj<-  Krig(x, Y, cov.function = stationary.cov, m = m, scale.type = scale.type, 
                method=method, GCV = GCV,
                cov.args = list(Covariance = "RadialBasis", 
                M = m, dimension = 2, Distance = "rdist.earth", 
                Dist.args = list(miles = miles)), ...)
              
    }
    obj$call<- match.call()
    class( obj) <- c("Krig", "Tps")
    return(obj)
}
