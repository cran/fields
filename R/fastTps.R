# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"fastTps" <-
function (x, Y, m = NULL, p = NULL,  theta,
    ...) 
{

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
    Tpscall$cov.function <- 
            "Compactly supported Wendland function "
    mKrig(x, Y, cov.function = "wendland.cov", 
          m = m, cov.args=list(k=p),  ... )
}

