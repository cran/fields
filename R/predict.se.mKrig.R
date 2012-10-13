# fields, Tools for spatial data
# Copyright 2004-2009, Institute for Mathematics Applied to Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"predict.se.mKrig" <- function(object, x = NULL, verbose = FALSE, 
    ...) {
    #
    # name of covariance function
    call.name <- object$cov.function.name
    #
    # default is to predict at data x's
    if (is.null(x)) {
        x <- object$x
    }
    x <- as.matrix(x)
    if (verbose) {
        print(x)
    }
    lambda <- object$lambda
    rho <- object$rhohat
    sigma2 <- lambda * rho
    if (verbose) {
        print(c(lambda, rho, sigma2))
    }
    k0 <- do.call(call.name, c(object$args, list(x1 = object$x, 
        x2 = x)))
    t0 <- t(fields.mkpoly(x, m = object$m))
    #
    # old form based on the predict function
    #   temp1 <-  rho*(t0%*% object$Omega %*%t(t0)) -
    #          rho*predict( object, y= k0, x=x) -
    #          rho*predict( object, y= k0, x=x, just.fixed=TRUE)

    # alternative formula using the d and c coefficients directly. 
    hold <- mKrig.coef(object, y = k0)
    temp1 <- rho * (colSums(t0 * (object$Omega %*% t0)) - colSums((k0) * 
        hold$c) - 2 * colSums(t0 * hold$d))
    # find marginal variances -- trival in the stationary case!
    temp0 <- rho * do.call(call.name, c(object$args, list(x1 = x, 
        marginal = TRUE)))
    # Add marginal variance to part from estimate
    temp <- temp0 + temp1
    return(sqrt(temp))
}
