# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"setup.image.smooth" <- function(nrow = 64, ncol = 64, 
    dx = 1, dy = 1, kernel.function = double.exp, theta = 1, 
    xwidth = nrow * dx, ywidth = ncol * dx, lambda = NULL, ...) {
    M2 <- round((nrow + xwidth/dx)/2)
    N2 <- round((ncol + ywidth/dy)/2)
    M <- 2 * M2
    N <- 2 * N2
    xi <- seq(-(M2 - 1), M2, 1) * dx
    xi <- xi/theta
    
    yi <- seq(-(N2 - 1), (N2), 1) * dy
    yi <- yi/theta
    dd <- sqrt((matrix(xi, M, N)^2 + matrix(yi, M, N, byrow = TRUE)^2))
    out <- matrix(kernel.function(dd, ...), nrow = M, ncol = N)
    out2 <- matrix(0, M, N)
    out2[M2, N2] <- 1
    
    W = fft(out)/fft(out2)
    if (!is.null(lambda)) {
        # want fft(out) / ( fft(out2)*lambda + fft(out))
        W = W/(lambda/fft(out2) + W)
    }
    
    list(W = W/(M * N), dx = dx, dy = dy, xwidth = xwidth, ywidth = ywidth, 
        M = M, N = N, m = nrow, n = ncol, lambda = lambda, grid = list(x = xi, 
            y = yi))
}
