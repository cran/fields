"matern.cov" <-
function (x1, x2, theta = rep(1, ncol(x1)), smoothness = 0.5) 
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
    matrix(matern(c(rdist(x1, x2)), smoothness = smoothness)$y, 
        nrow = n1, ncol = n2)
}
