# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
matern.image.cov <- function(ind1, ind2, Y, cov.obj = NULL, 
    setup = FALSE, grid, M = NULL, N = NULL,theta=1, smoothness=.5 ) {
    if (is.null(cov.obj)) {
        dx <- grid$x[2] - grid$x[1]
        dy <- grid$y[2] - grid$y[1]
        m <- length(grid$x)
        n <- length(grid$y)
        if (is.null(M)) 
            M <- ceiling2(2 * m)
        if (is.null(N)) 
            N <- ceiling2(2 * n)
# make sure M and N are even.
# (not sure what it means if this is not the case!)        
        if( M%%2 !=0) {
              M<- M+1}
        if( N%%2 !=0) {
              N<- N+1}    
# need to evaluate the covariance between the center of the grid and
# every grid point do this using several simple steps for efficiency.
        xGrid<- (1:M) * dx - (dx * M)/2
        yGrid<- (1:N) * dy -  (dy * N)/2
# a matrix the same size as the grid that has the distance between every
# grid point and the center point. 
        bigDistance<-
           sqrt(
             matrix( xGrid^2, M,N, byrow=FALSE) +
             matrix( yGrid^2, M,N, byrow=TRUE) )
# this should make for a nice image plot of the covariance w/r to the center point #       
        out<- Matern( bigDistance /theta, smoothness=smoothness)
        temp <- matrix(0, nrow = M, ncol = N)
        temp[M/2, N/2] <- 1
        wght <- fft(out)/(fft(temp) * M * N)
        cov.obj <- list(m = m, n = n, grid = grid, N = N, M = M, 
            wght = wght, call = match.call())
        if (setup) {
            return(cov.obj)
        }
    }
    temp <- matrix(0, nrow = cov.obj$M, ncol = cov.obj$N)
    if (missing(ind1)) {
        temp[1:cov.obj$m, 1:cov.obj$n] <- Y
        Re(fft(fft(temp) * cov.obj$wght, inverse = TRUE)[1:cov.obj$m, 
            1:cov.obj$n])
    }
    else {
        if (missing(ind2)) {
            temp[ind1] <- Y
        }
        else {
            temp[ind2] <- Y
        }
        Re(fft(fft(temp) * cov.obj$wght, inverse = TRUE)[ind1])
    }
}
