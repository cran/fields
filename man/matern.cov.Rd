\name{matern.cov}
\alias{matern.cov}
\title{
  Matern covariance function 
}
\description{
Given two sets of locations computes the Matern cross covariance matrix
for covariances among all pairings. 
}
\usage{
matern.cov
(x1, x2, theta = rep(1, ncol(x1)), smoothness = 0.5)
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
Range (or scale) parameter. This can be a scalar or a vector that is the
same length as the dimension of the locations. Default is theta=1.
}
\item{smoothness}{
The shape parameter for the Matern family. The exponential
is found with smoothness = 0.5 as smoothness goes to infinity one
recovers the Gaussian. 
}
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
In d dimensions a process with Matern smoothness parameter S
will have S + d/2 derivatives that exist in a mean square sense. 

Functional Form: If x1 and x2 are matrices where nrow(x1)=m and nrow(
x2)=n then this
function should return a mXn matrix where the (i,j) element is the
covariance between the locations x1[i,] and x2[j,]. The
covariance is found as  H( D.ij) where  D.ij is the Euclidean
distance between  x1[i,] and x2[j,] but having first been scaled by theta.
H is proportional to a modified Bessel function that depends on the
smoothness. 
H is normalized so that H(0)=1.  
Specifically the definition of the distance matrix is 

D.ij = sqrt(  sum.k (( x1[i,k] - x2[j,k]) /theta[k])**2 ).

Note that if theta is a scalar then this defines an isotropic covariance
function. 

Implementation: The function rdist is a useful FIELDS function that finds
the cross
distance matrix ( D defined above) for two sets of locations. Thus in
compact S code we have  

u <- t(t(x1)/theta)

v <- t(t(x2)/theta)

H(-rdist(u, v))

FORTRAN: The actual function calls FORTRAN to evaluate the special
function H (the FIELDS function matern). We
use the code by Montse Fuentes for this purpose. 
 
A simple modification of this function for the user would be to substitute
rdist.earth for rdist to give a distance metric that makes sense for
lon/lat coordinates. 
}
\seealso{
 Krig, matern.cov, rdist, rdist.earth, gauss.cov, exp.image.cov 
}
\examples{
# maternl covariance matrix ( marginal variance =1) for the ozone
#locations 
out<- matern.cov( ozone$x, theta=100)

# out is a 20X20 matrix
out2<- matern.cov( ozone$x[6:20,],ozone$x[1:2,], theta=100)
# out2 is 15X2 matrix 
# Kriging fit where the nugget variance is found by GCV 
fit<- Krig( ozone$x, ozone$y, matern.cov, theta=100, smoothness=.8)

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
