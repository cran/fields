"plot.sim.krig.image" <-
function (obj) 
{
    par.old <- par()
    on.exit(par(old.par))
    M <- length(obj$out)
    n <- round(sqrt(M))
    m <- round(M/n)
    if (m * n < M) 
        n <- n + 1
    set.panel(m, n, relax = T)
    for (k in 1:M) {
        image.plot(x = obj$grid$x, y = obj$grid$y, z = obj$out[[k]])
    }
}
