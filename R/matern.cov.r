"matern.cov" <-
function (x1, x2, theta = rep(1, ncol(x1)), smoothness = 0.5, 
    scale = 1) 
{
    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if (missing(x2)) 
        x2 <- x1
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)

# convert single range parameters or just diagonals
# into a matrix
    if (length(theta) == 1)         theta <- rep(theta, ncol(x1))
    if( is.vector(theta)) theta<- diag( theta)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
#OLD    x1 <- t(t(x1)/theta)
#OLD    x2 <- t(t(x2)/theta)
     
    x1 <- x1%*% t(solve( theta))
    x2 <- x2%*% t(solve( theta))

    matrix(matern(c(rdist(x1, x2)), smoothness = smoothness, 
        scale = scale), nrow = n1, ncol = n2)
}
