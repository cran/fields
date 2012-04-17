# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
cubic.cov <- function(x1, x2, theta = 1, C = NA, marginal = FALSE) {
    # comments in Exp.simple.cov for more details about the
    # required parts of this covariance
    
    if (is.matrix(x1)) {
        if (ncol(x1) != 1) {
            stop(" x is a matrix this is a  1-d covariance")
        }
    }
    
    if (missing(x2)) {
        x2 <- x1
    }
    # local function
    fun.temp <- function(u, v) {
        1 + ifelse(u < v, v * (u^2)/2 - (u^3)/6, u * (v^2)/2 - 
            (v^3)/6)
    }
    if (is.na(C[1]) & !marginal) {
        # cross covariance matrix
        return(outer(c(x1), c(x2), FUN = fun.temp))
    }
    if (!is.na(C[1])) {
        # product of cross covariance with a vector
        return(outer(c(x1), c(x2), FUN = fun.temp) %*% C)
    }
    if (marginal) {
        # marginal variance
        return((x1^3)/3)
    }
}
