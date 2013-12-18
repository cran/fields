# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"stats" <- function(x, by) {
    if (!missing(by)) {
        x <- cat.to.list(c(x), by)
    }
    if (!is.list(x) & !is.matrix(x)) 
        x <- matrix(x, ncol = 1)
    if (is.list(x)) {
        ncol <- length(x)
        out <- matrix(NA, ncol = ncol, nrow = length(describe()))
        dimnames(out) <- list(describe(), names(x))
        for (j in (1:ncol)) {
            if (is.numeric(x[[j]])) {
                out[, j] <- describe(x[[j]])
            }
        }
        return(out)
    }
    if (is.matrix(x)) {
        nc <- ncol(x)
        out <- matrix(NA, ncol = nc, nrow = length(describe()))
        dimnames(out) <- list(describe(), dimnames(x)[[2]])
        for (j in (1:nc)) {
            out[, j] <- describe(x[, j])
        }
        return(out)
    }
}
