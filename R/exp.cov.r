"exp.cov" <-
function (x1, x2, theta = rep(1, ncol(x1)), p = 1, C = NA) 
{
    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if (missing(x2)) 
        x2 <- x1
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)
    if (length(theta) == 1) 
        theta <- rep(theta, ncol(x1))
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    x1 <- t(t(x1)/theta)
    x2 <- t(t(x2)/theta)
    par <- p
    if (is.na(C[1])) {
        exp(-rdist(x1, x2)^p)
    }
    else {
        .Fortran("multeb", nd = as.integer(d), x1 = as.double(x1), 
            n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
            par = as.double(par), c = as.double(C), h = as.double(rep(0, 
                n1)), work = as.double(rep(0, n2)),PACKAGE="fields")$h
    }
}
