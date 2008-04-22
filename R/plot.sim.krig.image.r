# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"plot.sim.krig.image" <-
function (x,...) 
{
    obj<- x # hack S3
    par.old <- par(read.only=TRUE)
    on.exit(par(par.old))
    M <- length(obj$out)
    n <- round(sqrt(M))
    m <- round(M/n)
    if (m * n < M) 
        n <- n + 1
    set.panel(m, n, relax = TRUE)
    for (k in 1:M) {
        image.plot(x = obj$grid$x, y = obj$grid$y, z = obj$out[[k]])
    }
}
