# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"Rad.cov" <- function(x1, x2, p = 1, m=NA,with.log = TRUE, 
    with.constant = TRUE, C = NA, marginal = FALSE, derivative = 0) {
    #
    # mth order thin plate spline radial basis functions
    # in d dimensions
    # usually called with p 2m-d
    #  marginal dummy argument
    #  this should only be called within predict.se.Krig
    #  and provides the correct calculation. Because this is
    #  a generalized covariance the marginal variance is not really
    #  defined.
    #
    if (marginal) {
        return(rep(0, nrow(x1)))
    }
    #
    # coerce locations to matrices, if x2 is missing use x1
    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    if( is.na(m)){
      m <- (d + p)/2}
    else{
      p <- 2*m -d}
    if( p < 0 ) {
      stop(" p is negative (m possibly too small)")}
    # parameter list to send to the FORTRAN
      par <- c(p/2, ifelse((d%%2 == 0) & (with.log), 1, 0))
    #
    # multiply by constant if requested
      rbf.constant <- ifelse(with.constant, radbas.constant(m, 
        d), 1)
    # compute matrix in FORTRAN
    if (is.na(C[1])) {
        temp <- .Fortran("radbas", nd = as.integer(d), x1 = as.double(x1), 
            n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
            par = as.double(par), k = as.double(rep(0, n1 * n2)))
        return(rbf.constant * matrix(temp$k, ncol = n2, nrow = n1))
    }
    else {
        #   do cross covariance matrix multiplication in FORTRAN
        if (derivative == 0) {
            #     evaluate function not partial derivatives.
            C <- as.matrix(C)
            n3 <- ncol(C)
            temp <- .Fortran("multrb", nd = as.integer(d), x1 = as.double(x1), 
                n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
                par = as.double(par), c = as.double(C), n3 = as.integer(n3), 
                h = as.double(rep(0, n1 * n3)), work = as.double(rep(0, 
                  n2)))$h
            return(rbf.constant * matrix(temp, nrow = n1, ncol = n3))
        }
        else {
            if (ncol(C) > 1) {
                stop("Can only evaluate derivatives on one spline fit")
            }
            temp <- .Fortran("mltdrb", nd = as.integer(d), x1 = as.double(x1), 
                n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
                par = as.double(par), c = as.double(C), h = as.double(rep(0, 
                  n1 * d)), work = as.double(rep(0, n2)))$h
            return(rbf.constant * matrix(temp, nrow = n1, ncol = d))
        }
    }
    stop("should not get here!")
}
