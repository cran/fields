"sim.rf.W" <-
function (cov.obj) 
{
    N <- cov.obj$n
    M <- cov.obj$m
    z <- sqrt(cov.obj$D) * matrix(rnorm(N * M), ncol = N, nrow = M)
    Wtransform.image(z, inv = T, cut.min = cov.obj$cut.min)
}
