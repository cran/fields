# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"setup.image.smooth" <-
function ( nrow = 64, ncol = 64, dx = 1, dy = 1, kernel.function = 
double.exp, 
    theta = 1, xwidth = nrow*dx, ywidth = ncol*dx, ...) 
{
    M2 <- round((nrow + xwidth/dx)/2)
    N2 <- round((ncol + ywidth/dy)/2)
    M <- 2 * M2
    N <- 2 * N2
    xi <- (0:(M2-1)) * dx
    xi <- xi/theta
    yi <- (0:(N2-1)) * dy
    yi <- yi/theta
    out <- kernel.function(
             (matrix(xi, M2, N2)^2 + matrix(yi,M2, N2, byrow = TRUE)^2), ...)/theta
    out <- cbind(out, out[, N2:1])
    out <- rbind(out, out[M2:1, ])
    
    list( W=fft(out)/(M * N), dx=dx, dy=dy,xwidth=xwidth, ywidth=ywidth)
}
