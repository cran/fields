"matern.cov" <-
function (x1, x2, theta = 1.0, smoothness = 0.5, 
    scale = 1) 
{
# coerce to matrix
    if (!is.matrix(x1)) 
        x1 <- matrix(c(x1), ncol=1)

    if (missing(x2)) 
        x2 <- x1
    if (!is.matrix(x2)) 
        x2 <- matrix(c(x2), ncol=1)

    if (length(theta) == 1) 
        theta <- rep(theta, ncol(x1))

# handle special case of 1-d
    if( ncol(x1)==1) { theta<- matrix( c(theta),1,1)}

# handle special case of just diagonal elements of  theta
    if (is.vector(theta)) 
        theta <- diag(theta)

# following now treats theta as a full matrix for scaling and rotation. 

    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    x1 <- x1 %*% t(solve(theta))
    x2 <- x2 %*% t(solve(theta))

    matrix(matern(c(rdist(x1, x2)), smoothness = smoothness, 
        scale = scale), nrow = n1, ncol = n2)
}

