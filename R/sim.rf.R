# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"sim.rf" <-
function (obj) 
{
    n <- obj$n
    m <- obj$m
    M <- obj$M
    N <- obj$N
    if (any(Re(obj$wght) < 0)) {
        stop("FFT of covariance has negative\nvalues")
    }
    z <- fft(matrix(rnorm(N * M), ncol = N, nrow = M))
    Re(fft(sqrt(obj$wght) * z, inverse = TRUE))[1:m, 1:n]/sqrt(M * 
        N)
}
