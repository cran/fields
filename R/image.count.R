"image.count" <-
function (x, grid = NULL, nrow = 64, ncol = 64) 
{
    Z <- rep(1, nrow(x))
    if (is.null(x) & is.null(grid)) {
        grid <- list(x = 1:nrow, y = 1:ncol)
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
    rep.info <- cat.matrix(ind)
    uniquerows <- !dup(rep.info)
    if (sum(uniquerows) < length(Z)) {
        ind <- ind[uniquerows, ]
        Z <- fast.1way(rep.info, Z)$n
    }
    m <- length(grid$x)
    n <- length(grid$y)
    temp <- matrix(NA, nrow = m, ncol = n)
    temp[ind] <- Z
    call <- match.call()
    list(x = grid$x, y = grid$y, z = temp, call = call)
}
