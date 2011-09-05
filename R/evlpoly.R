# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
fields.evlpoly <- function(x, coef) {
    # evaluates polynomial at x values with coefficients coef[i] and powers i-1
    #
    n <- length(x)
    J <- length(coef)
    results <- rep(0, n)
    temp <- .Fortran("evlpoly", x = as.double(x), n = as.integer(n), 
        coef = as.double(coef), j = as.integer(J), results = as.double(results) 
       )$results
    return(temp)
}
