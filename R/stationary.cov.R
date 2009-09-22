# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"stationary.cov" <- function(x1, x2, Covariance = "Exponential", 
    Distance = "rdist", Dist.args = NULL, theta = 1, V = NULL, 
    C = NA, marginal = FALSE, ...) {
    # get covariance function arguments from call
    Cov.args <- list(...)
    # coerce x1 and x2 to matrices
    if (is.data.frame(x1)) 
        x1 <- as.matrix(x1)
    if (!is.matrix(x1)) 
        x1 <- matrix(c(x1), ncol = 1)
    if (missing(x2)) 
        x2 <- x1
    if (is.data.frame(x2)) 
        x2 <- as.matrix(x2)
    if (!is.matrix(x2)) 
        x2 <- matrix(c(x2), ncol = 1)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    #
    # separate out a single scalar transformation and a
    # more complicated scaling and rotation.
    # this is done partly to have the use of great circle distance make sense
    # by applying the scaling  _after_ finding the distance.
    #
    # tranform coordinates if theta not a scalar
    if (length(theta) > 1) {
        stop("theta as a vector matrix has been depreciated use the V argument")
    }
    #
    # following now treats V as a full matrix for scaling and rotation.
    #
    # try to catch incorrect conbination  of great circle distance and V
    if (Distance == "rdist.earth" & !is.null(V)) {
        stop("V not supported with great circle distance")
    }
    if (!is.null(V)) {
        if (theta != 1) {
            stop("can't specify both theta and V!")
        }
        x1 <- x1 %*% t(solve(V))
        x2 <- x2 %*% t(solve(V))
    }
    #
    # locations are now scaled and rotated correctly
    # now apply covariance function to pairwise distance matrix, or multiply
    # by C vector or just find marginal variance
    # this if block finds the cross covariance matrix
    if (is.na(C[1]) & !marginal) {
        # note overall scalling by theta (which is just theta under isotropic case)
        bigD <- do.call(Distance, c(list(x1 = x1, x2 = x2), Dist.args))/theta
        return(do.call(Covariance, c(list(d = bigD), Cov.args)))
    }
    # or multiply cross covariance by C
    if (!is.na(C[1])) {
        #
        # as coded below this is not particularly efficient of memory
        #
        bigD <- do.call(Distance, c(list(x1 = x1, x2 = x2), Dist.args))/theta
        return(do.call(Covariance, c(list(d = bigD), Cov.args)) %*% 
            C)
    }
    # or find marginal variance and return  a vector.
    if (marginal) {
        sigma2 <- do.call(Covariance, c(list(d = 0), Cov.args))
        return(rep(sigma2, nrow(x1)))
    }
    # should not get here
}
