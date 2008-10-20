# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"image.smooth" <-
function (x, wght = NULL, dx = 1, dy = 1,         
       kernel.function = double.exp, theta = 1, grid = NULL, tol = 1e-08,
       xwidth = NULL, ywidth = NULL, 
       weights = NULL,...) 
{
# first part of this function is figuring what has been passed and
# what to do
  
    if( is.list( x)) {
# assume that an image list format has been passed as x
      Y<- x$z
     grid<- list( x= x$x, y=x$y)}
    else{   
     Y<- x} # hack S3
    
    if (!is.matrix(Y)) {
        stop("Requires a matrix")
    }

    m <- nrow(Y)
    n <- ncol(Y)

# use information in previous setup kernel function from a
# a call to setup.image.smooth and in the process override any
# passed arguments
    
    if( !is.null(wght)){
      dx<- wght$dx
      dy<- wght$dy
      xwidth<- wght$xwidth
      ywidth<- wght$ywidth}
      
# set up grid if it is missing    
    if (is.null(grid)) {
     grid<- list( x= (1:m)*dx, y=(1:n)*dy)}
    else{
        dx <- grid$x[2] - grid$x[1]
        dy <- grid$y[2] - grid$y[1]}

# padding of zeroes around actual image
# if less than m and n there may be spurious effects due to
# the periodicity from the fft.
# make sure that the span of the kernel is less than xwidth and ywidth.
# there will be substantial speedup if the kernel has a small support,     
# Y is big  (e.g. 512X512) and Mwidth and N widht are adjusted to suit.
    
    if (is.null(xwidth)){
       xwidth<- dx*m}
  
    if (is.null(ywidth)){
       ywidth<- dy*n}
  
# kernel wght function as fft
# reusing this saves an fft for each image smooth.     
    if (is.null(wght)) {
      wght <- setup.image.smooth(nrow = m, ncol = n, xwidth = xwidth, 
              ywidth = ywidth, dx = dx, dy = dy,
              kernel.function = kernel.function, theta = theta)
    }

    M <- nrow(wght$W)
    N <- ncol(wght$W)
    temp <- matrix(0, nrow = M, ncol = N)
    temp2 <- matrix(0, nrow = M, ncol = N)
# pad with zeroes     
    if (!is.null(weights)) {
        temp[1:m, 1:n] <- Y * weights
        temp[is.na(temp)] <- 0
        temp2[1:m, 1:n] <- ifelse(!is.na(Y), weights, 0)
    }
    else {
        temp[1:m, 1:n] <- Y
        temp[is.na(temp)] <- 0
        temp2[1:m, 1:n] <- ifelse(!is.na(Y), 1, 0)
    }
# temp and temp2 are numerator and denominator of Nadarya-Watson estimator.    
    temp <- Re(fft(fft(temp) * wght$W, inverse = TRUE))[1:m, 1:n]
    temp2 <- Re(fft(fft(temp2) * wght$W, inverse = TRUE))[1:m,1:n]
    
# try not to divide by zero!    
    temp <- ifelse((temp2 > tol), (temp/temp2), NA)
       
    list(x = grid$x, y = grid$y, z = temp)

}
