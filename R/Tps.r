# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"Tps" <- function(x, Y, m = NULL, p = NULL, scale.type = "range", 
    ...) {
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
    Tpscall <- match.call()
    Tpscall$cov.function <- "Thin plate spline radial basis functions (Rad.cov) "
    Krig(x, Y, cov.function = Rad.cov, m = m, scale.type = scale.type, 
        outputcall = Tpscall, p = p, GCV = TRUE, ...)
}
