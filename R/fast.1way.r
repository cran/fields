"fast.1way" <-
function (lev, y, w = rep(1, length(y))) 
{
    N <- length(y)
    tags <- lev[!dup(lev)]
    lev <- match(lev, tags)
    id <- order(lev)
    brk <- c(diff(lev[id]) != 0, TRUE)
    w.means <- diff(c(0, cumsum(w[id])[brk]))
    means <- diff(c(0, cumsum(y[id] * w[id])[brk]))/w.means
    n <- diff(c(0, (1:N)[brk]))
    SSE <- sum(w * (y - means[lev])^2)
    MSE <- SSE/(length(y) - length(n))
    names(means) <- tags
    names(w.means) <- tags
    list(means = means, SSE = SSE, w.means = w.means, n = n, 
        MSE = MSE, fitted.values = means[lev], tags = tags)
}
