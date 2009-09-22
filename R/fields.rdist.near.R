# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
fields.rdist.near <- function(x1, x2, delta, max.points = NULL, 
    mean.neighbor = 50) {
    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if (missing(x2)) 
        x2 <- x1
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    if (is.null(max.points)) {
        Nmax <- n1 * mean.neighbor
    }
    else {
        Nmax <- max.points
    }
    out <- .Fortran("ddfind", nd = as.integer(d), x1 = as.double(x1), 
        n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
        D0 = as.double(delta), ind = as.integer(rep(0, Nmax * 
            2)), rd = as.double(rep(-1, Nmax)), Nmax = as.integer(Nmax), 
        iflag = as.integer(1), PACKAGE = "fields")
    N <- out$Nmax
    if (out$iflag == -1) {
        cat("temp space set at", Nmax, fill = TRUE)
        stop("Ran out of space, increase max.points")
    }
    return(list(ind = matrix(out$ind, Nmax, 2)[1:N, ], ra = out$rd[1:N], 
        da = c(n1, n2)))
}
