\name{interp.surface}
\alias{interp.surface}
\title{
  Fast bilinear interpolator from a grid.
}
\description{
Uses bilinear weights to interpolate values on a rectangular
grid to arbitrary locations. 
}
\usage{
interp.surface(obj, loc)
}
\arguments{
\item{obj}{
A list with components x,y, and z in the same style as used by contour,
persp, image etc. x and y are the X and Y grid values and z is a matrix
with the corresponding values of the surface
}
\item{loc}{
A matrix of (irregular) locations to interpolate. First column of loc
isthe X coordinates and second is the Y's. 
}
}
\value{
 An vector of interpolated values.  NA are returned for
regions of the obj\$z that are NA and also for locations outside of the
range of the grid. 
}
\details{
Here is a brief explanation of
the interpolation:  Suppose that the location, (locx, locy)  lies in
between the first two grid points in both x an y. That is  locx is between
x1 and x2 and
locy is between y1 and y2.  Let ex= (l1-x1)/(x2-x1) ey= (l2-y1)/(y2-y1).
The
interpolant is

( 1-ex)(1-ey)*z11 + (1- ex)(ey)*z12 + ( ex)(1-ey)*z21 + ( ex)(ey)*z22  

Where the z's are the corresponding elements of the Z matrix. 
}
\seealso{
image.smooth, as.surface, as.image, image.plot, krig.image,Tps
}
\examples{
#
# evaluate an image at a finer grid
# 

data( lennon)
# create the surface object
obj<- list( x= 1:20, y=1:20, z= lennon[ 201:220, 201:220])

# sample at 50 equally spaced points 
temp<- seq( 1,20,,50)
make.surface.grid( list( temp,temp))-> loc
interp.surface( obj, loc)-> look
# take a look
image.plot( as.surface( loc, look))

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
