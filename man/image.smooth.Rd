\name{image.smooth}
\alias{image.smooth}
\title{
  Kernel smoother for irregular 2-d data 
}
\description{
Takes an image matrix and applies a kernel smoother to it. Missing values
are handled using the Nadaraya/Watson normalization of the kernel. 
}
\usage{
image.smooth
 (Y, wght = NULL, dx = 1, dy = 1, Nwidth = nrow(Y), Mwidth = ncol(Y),
    kernel.function = function(x) {
        exp(-abs(x))
    }, theta = 1, grid = NULL, tol = 1e-08, xwidth = NULL, ywidth = NULL,
    weights = NULL)
}
\arguments{
\item{Y}{
A matrix image with missing values are indicated by NA. 
}
\item{wght}{
FFT of smoothing kernel. If this is NULL the default is to compute this
object. 
}
\item{grid}{
A list with x and y components. Each are equally spaced and define the rectangular. ( see grid.list)
}
\item{dx}{
Grid spacing in x direction 
}
\item{dy}{
Grid spacing in x direction 
}
\item{Nwidth}{
Half the width of the kernel in the y direction in numbers of grid 
points.  
Default is half the number of columns of Y 
}

\item{Mwidth}{
Half the width of the kernel in the x direction in numbers of grid 
points.
Default is half the number of columns of Y 
}
\item{kernel.function}{
An S function that takes as its argument the squared distance between two 
points divided by the bandwidth. The default is exp( -x) yielding a normal 
kernel 
}
\item{theta}{
The bandwidth 
}
\item{\dots}{
Other arguments to be passed to the kernel function 
}
\item{xwidth}{
This has the same function as Mwidth but is in the scale of the x grid.
If specified  then  Mwidth = xwidth/ ( grid$x[2]- grid$x[1])
}
\item{ywidth}{
This has the same function as Nwidth but is in the scale of the x grid.
If specified  then  Nwidth = xwidth/ ( grid$x[2]- grid$x[1])
}
\item{weights}{
Weights to apply when smoothing.}
\item{tol} {
Tolerance for the weights of the N-W kernel. This avoids kernel estimates that are "far" away from data.
Grid points with weights less than tol are set to NA.
 }


}
\value{
The smoothed matrix. 
}
\details{
The function works by taking convolutions using an FFT. The missing pixels 
are taken into account and the kernel smoothing is correctly normalized
following the classical Nadaraya-Watson estimator. If the kernel has 
limited support then the width arguments can be set to reduce the amount 
of computation. (See example below.) }
\seealso{
image.smooth.setup, as.image  
}
\examples{
# first convert precip data to the 128X128 discretized image format ( with 
# missing  values to indicate where data is not observed) 
# 
out<- as.image( precip$y, x= precip$x, nrow=128, ncol=128) 
# out$z is the image matrix 

dx<- out$x[2]- out$x[1] 
dy<-  out$y[2] - out$y[1] 
#  
look<- image.smooth( out$z, dx=dx, dy=dy, theta= .25) 

# grid scale in degree and so kernel bandwidth is .25 degrees. 
image.plot( x= out$x, y=out$y, z= look) 
points( precip$x)

# to save on computation, decrease the padding with zeroes 

look<- image.smooth( out$z, dx=dx, dy=dy, theta= .25, Mwidth=32,Nwidth=32) 

# the range of these data is ~ 10 and so  32*( 10/128) =  2.5 
# about 10 standard deviations of the normal kernel so there is still 
# lots of room for padding  
# a minimal choice might be  Mwidth = 4* (.25/dx)  4 SD for the normal  
# 
# creating weighting object outside the call  
# this is useful when one wants to smooth different data sets but on the 
# same grid with the same kernel function 
# 

wght<- image.smooth.setup( nrow=128, ncol=128,  dx=dx, dy=dy, theta= .25,
Mwidth=32, Nwidth=32)

#
# 6 random fields from smoothing white noise with this filter.
#
set.panel( 3,2)
par( pty="s")
look<- image.smooth( out$z, wght) 
for( k in 1:6) {
look<- image.smooth( matrix( rnorm(128**2), 128,128), wght)
image( look)
}
}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
