"setup.W.image.cov" <-
function (m, n, weights = 1, D = NULL, cut.min = 8, max.m = m, 
    max.n = n, x = NULL, grid = NULL, normalize = T, ...) 
{
    call <- match.call()
    if (!is.null(x)) {
        grid <- list(x = 1:m, y = 1:n)
    }
    if (!is.null(grid)) {
        m <- length(grid$x)
        n <- length(grid$y)
    }
    if (is.null(D)) {
        junk <- Wtransform.D(m, n, weights = weights, cut.min = cut.min)
    }
    else {
        junk <- list(D = D, max.m = max.m, max.n = max.n)
    }
    e <- matrix(0, m, n)
    x1 <- m/2
    x2 <- n/2
    if (normalize) {
        e[x1, x2] <- 1
        temp.max <- max(Wtransform.image(junk$D * Wtransform.image(e, 
            inv = T, transpose = T, cut.min = cut.min), inv = T, 
            cut.min = cut.min, ))
    }
    else {
        temp.max <- 1
    }
    return(list(m = m, n = n, D = junk$D/temp.max, grid = grid, 
        cut.min = cut.min, max.n = junk$max.n, max.m = junk$max.m, 
        weights = weights, call = call))
}
