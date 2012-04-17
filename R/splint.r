# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"splint" <- function(x, y, xgrid, wt = NULL, derivative = 0, 
    lam = 0, df = NA, lambda = NULL) {
    #
    # reform calling args if passed as a matrix or list
    
    if (is.matrix(x)) {
        if (ncol(x) > 1) {
            xgrid <- y
            y <- x[, 2]
            x <- x[, 1]
        }
    }
    if (is.list(x)) {
        xgrid <- y
        y <- x$y
        x <- x$x
    }
    if (any(duplicated(x))) {
        stop("duplicated x values, use sreg")
    }
    if ((derivative > 2) | (derivative < 0)) 
        stop("derivative must be 0,1,2")
    if (length(x) != length(y)) 
        stop("Lengths of x and y must match")
    n <- length(x)
    #default values for weights
    # NOTE: weights do not matter when interpolating (lam==0)
    if (is.null(wt)) {
        wt <- rep(1, n)
    }
    # find lambda from eff degrees of freedom if it is passed
    if (!is.na(df)) {
        if ((df < 2) | (df > n)) {
            stop("df out of range")
        }
        lam <- sreg.df.to.lambda(df, x, wt)
    }
    # use lambda is it is passed
    if (!is.null(lambda)) {
        lam <- lambda
    }
    igcv <- ifelse(lam == 0, 2, 0)
    # call to old FORTRAN
    temp <- .Fortran("css", h = as.double(ifelse(igcv == 2, 1, 
        log(lam))), as.integer(n), as.double(x), as.double(y), 
        wt = as.double(1/sqrt(wt)), sy = as.double(rep(0, n)), 
        as.double(1), as.double(1), as.double(1), as.integer(length(xgrid)), 
        as.double(xgrid), ygrid = as.double(rep(0, length(xgrid))), 
        job = as.integer(c(igcv, 3, 0)), as.integer(derivative), 
        as.integer(0))
    return(temp$ygrid)
}
