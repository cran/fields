\name{matern.cov}
\alias{matern.cov}
\alias{matern.earth.cov}
\title{
  Matern covariance function 
}
\description{
Given two sets of locations computes the Matern cross covariance matrix
for covariances among all pairings. In the case of matern.earth.cov the 
locations are assumed to be lon/lat coordinates and  the covariance 
function is a spherical modification. 
 }
\usage{
matern.cov(x1, x2, theta = 1.0, smoothness = 0.5, scale=1)
matern.earth.cov(x1, x2, theta = 1, 
smoothness = 0.5, scale = 1, miles = TRUE, R = NULL)
}
\arguments{
\item{x1}{
Matrix of first set of locations where each row gives the coordinates of a
particular
point.
}
\item{x2}{
Matrix of second set of locations where each row gives the coordinates of
a particular point. If this is missing x1 is used. 
}
\item{theta}{
Range (or scale) parameter matrix. This can be a scalar, a vector or 
matrix for matern.cov but must be a scalar for the spherical version.
 
Default is theta=1. The locations are scaled by the inverse of theta before the 
distances are found. Off diagonal elements of theta build in anisotrophy.
}
\item{smoothness}{
The shape parameter for the Matern family. The exponential
is found with smoothness = 0.5 as smoothness goes to infinity one
recovers the Gaussian. 
}
\item{scale}{ 
Marginal variance. 
}
  \item{miles}{If TRUE then radius of earth is in miles, if FALSE then 
radius is in kilometers. } 
  \item{R}{ Radius used for sphere default is to use the earth's 
radius. }

}
\value{
The cross covariance matrix between locations x1 and x2. If x1 is
equal to x2 then this is the covariance matrix for this set of locations. 
In general if nrow(x1)=m and nrow(
x2)=n then the returned matrix, Sigma will be mXn.
Under some cases when the range is very small  
calculation of the Bessel functions fails and if so an NA will be 
returned. 

}
\details{
Derivatives of sample paths:
In d dimensions a Gaussian process with Matern covariance and smoothness 
parameter 
nu will have nu + d/2 derivatives that exist in a mean square sense. 
The marginal variance is given by the parameter scale. When this 
parameter is set to one this yields a correlation 
function. 

There several different ways to parameterize the Matern family and the 
reader is referred to Stein's book page 49 for discussion. In terms 
of the more geostatistical terminology, we note that out theta is also the 
"range" and the scale is also the "sill" if there is no nugget variance 
included in the covariance. We caution that the range for this function 
gives a qualitatively different scaling as one varies the smoothness. 

Functional Form: If x1 and x2 are matrices where nrow(x1)=m and nrow(
x2)=n  and each row are the coordinates of a location, then this
function should return a mXn matrix where the (i,j) element is the
covariance between the locations x1[i,] and x2[j,]. The
covariance is found as  H( D.ij) where  D.ij is the Euclidean
distance between  x1[i,] and x2[j,] but having first been scaled by theta.
H is proportional to a modified Bessel function 
of third kind using denoted by K.nu .  In our parameterization we take 
smoothness = nu and H is normalized so that H(0)=1. (See the function matern for a succinct definition as R code.) 
The reader is referred 
to Stein's book, page 31 for more details. Note that we do not use Stein's 
normalization, however, as it seems more useful to normalize the 
covariance so that .     


Definition of the distance matrix: 
D.ij = sqrt(  sum.k (( x1[i,k] - x2[j,k]) /theta[k])**2 ).

Note that if theta is a scalar then this defines an isotropic covariance
function. 

Implementation: The function rdist is a useful FIELDS function that finds
the cross
distance matrix ( D defined above) for two sets of locations. Thus in
compact S code we have  

u <- t(solve(theta)%*%t(x1))
v <- t(solve(theta)%*%t(x2))
H(-rdist(u, v))

where \code{solve(theta)} is the (matrix) inverse for theta. 

The adaptation to the sphere is to transform the angular separation 
among locations by  2R sin( angle/2). This is a measure of distance  
is then used in the matern functional form for the covariance. 
For small regions this will approximate the usual Matern but is a valid 
covariance at all distances. (See last example below.)
}
\seealso{
 Krig, matern, rdist, rdist.earth, gauss.cov, exp.image.cov, matern.imag.cov 
}
\examples{
#
# Presenting the Matern family:
# the function matern is called by matern.cov
d<- seq( 0,5,,200)
sm<- seq( .5, 8,,5)
temp<- matrix( NA, 200, 5)
for ( k in 1:5){
temp[,k] <- matern(d, smoothness=sm[k])
}
matplot( d, temp, type="l", lty=1)
# note differing correlation scales depending on smoothness 


# Matern covariance matrix ( marginal variance =1) for the ozone
# locations 
out<- matern.cov( ozone$x, theta=100, smoothness=1.0)
# out is a 20X20 matrix

out2<- matern.cov( ozone$x[6:20,],ozone$x[1:2,], theta=100, 
smoothness=1.0)

# out2 is 15X2 cross covariance matrix 

# Kriging fit using a Matern covariance and where the nugget  and 
# sill variances are found by GCV 
fit<- Krig( ozone$x, ozone$y, matern.cov, theta=100, smoothness=1.0)


######## examples with  matern.earth.cov

# create a 1 degree patch of lon/lats  at the equator:
x<- seq( -1,1,,25)
y<- x
make.surface.grid( list( x=x, y=y))-> xg

# find covariance with the  (lon,lat) =(0,0)  
# theta = 30 , a degree at the equator is about 69.17333 miles. 

matern.earth.cov( xg,rbind( c(0,0)), theta = 30 )-> look

# covariance surface
image.plot( x,y,matrix(look, 25,25)) 

# compare to matern.cov( xg,rbind( c(0,0)), theta = 30/69.17333)

# now change the grid to be close to the north pole
# create a 1 degree patch of lon/lats  at  (0, 80N) :

x<-  seq( -1,1,,25)
y<- x + 80
make.surface.grid( list( x=x, y=y))-> xg

# evaluate the same covariance function at this point

matern.earth.cov( xg,rbind( c(0,80)), theta = 30 )-> look

# plot covariance surface 
#     stretching is an artifact of this Mercator projection 

image.plot( x,y,matrix(look, 25,25)) 

# a Kriging example with the ozone data

data(ozone2)

Krig( ozone2$lon.lat, ozone2$y[16,], cov.function= matern.earth.cov, 
    theta=300, na.rm=TRUE)-> look

surface( look, type="C")
US( add=TRUE, lwd=2)



}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
