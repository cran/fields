# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"qsreg.trace" <- function(x, y, lam, maxit = 50, maxit.cv = 10, 
    tol = 1e-04, offset = 0, sc = sqrt(var(y)) * 1e-07, alpha = 0.5, 
    wt = rep(1, length(x)), cost = 1) {
    N <- length(y)
    if (length(x) != length(y)) 
        stop(" X and Y do not match")
    h <- log(lam)
    temp <- .Fortran("rcss", h = as.double(log(lam)), npoint = as.integer(N), 
        x = as.double(x), y = as.double(y), wt = as.double(wt), 
        sy = as.double(rep(0, N)), trace = as.double(0), diag = as.double(rep(0, 
            N)), cv = as.double(0), ngrid = as.integer(0), xg = as.double(0), 
        yg = as.double(0), job = as.integer(c(3, 0, 0)), ideriv = as.integer(0), 
        din = as.double(c(cost, offset, maxit, tol, sc, alpha)), 
        dout = as.double(rep(0, 4)), ierr = as.integer(0))$dout
    return(temp[3])
}
