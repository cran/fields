"as.image" <-
function (Z, ind = NULL, grid = NULL, x = NULL, nrow = 64, ncol = 64, 
    weights = rep(1, length(Z)), na.rm=TRUE) 
{
    if (!is.null(ind)) 
        x <- ind
    if (is.null(x) & is.null(grid)) {
        grid <- list(x = 1:nrow, y = 1:ncol)
    }

if( na.rm){
ind2<- !is.na( Z)
Z<- Z[ind2] 
x<- x[ind2,]
weights<- weights[ind2]
}
    if (!is.null(x) & is.null(grid)) {
        temp <- Krig.discretize(x, nrow, ncol)
        grid <- temp$grid
        ind <- temp$index
    }
    if (!is.null(x) & !is.null(grid)) {
        temp <- Krig.discretize(x, grid = grid)
        ind <- temp$index
    }
    m <- length(grid$x)
    n <- length(grid$y)
    rep.info <- cat.matrix(ind)
    uniquerows <- !dup(rep.info)
    if (sum(uniquerows) < length(Z)) {
        ind <- ind[uniquerows, ]
        temp <- fast.1way(rep.info, Z, w = weights)
        Z <- temp$means
        Ncell <- temp$n
        temp2 <- matrix(0, nrow = m, ncol = n)
        temp2[ind] <- Ncell
        temp3 <- matrix(NA, nrow = m, ncol = n)
        temp3[ind] <- temp$w.means
    }
    else {
        temp2 <- matrix(0, nrow = m, ncol = n)
        temp2[ind] <- 1
        temp3 <- matrix(NA, nrow = m, ncol = n)
        temp3[ind] <- 1
    }
    temp <- matrix(NA, nrow = m, ncol = n)
    temp[ind] <- Z
    call <- match.call()
    list(x = grid$x, y = grid$y, z = temp, call = call, ind = ind, 
        N = temp2, weights = temp3)
}
