\name{exp.cov}
\alias{exp.cov}
\alias{exp.cov.simple}
\title{
  Exponential, Gaussian and "power" covariance family 
}
\description{
Given two sets of locations computes the cross covariance matrix for
covariances among all pairings. 
}
\usage{
exp.cov(x1, x2, theta = rep(1, ncol(x1)), p = 1, C = NA)
exp.cov.simple(x1, x2, theta =1)
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
\item{p}{
Exponent in the exponential form. p=1 gives an exponential and p=2 gives a
Gaussian. Default is the exponential form.
}
\item{C}{
A vector with the same length as the number of rows of x2. 
If specified the covariance matrix will be multiplied by this vector.
}
}
\value{
If the argument C is NULL the cross covariance matrix. Moreover if x1 is
equal to x2 then this is the covariance matrix for this set of locations. 
In general if nrow(x1)=m and nrow(
x2)=n then the returned matrix, Sigma will be mXn. 

If C is a vector of length n,  
then returned value is the multiplication of the cross covariance matrix
with this vector:  Sigma\%*\%C 
}
\details{
Functional Form: If x1 and x2 are matrices where nrow(x1)=m and nrow(
x2)=n then this
function should return a mXn matrix where the (i,j) element is the
covariance between the locations x1[i,] and x2[j,]. The
covariance is found as  exp( -(D.ij **p)) where  D.ij is the Euclidean
distance between  x1[i,] and x2[j,] but having first been scaled by theta.
Specifically 

D.ij = sqrt(  sum.k (( x1[i,k] - x2[j,k]) /theta[k])**2 ).

Note that if theta is a scalar then this defines an isotropic covariance
function. 

Implementation: The function r.dist is a useful FIELDS function that finds
the cross
distance matrix ( D defined above) for two sets of locations. Thus in
compact S code we have  

u <- t(t(x1)/theta)

v <- t(t(x2)/theta)

exp(-rdist(u, v)**p)

FORTRAN: The actual function calls FORTRAN to make the evaluation more
efficient this is especially important when the C argument is supplied.  
So unfortunately the actual code is not as crisp as the few lines given
above. For purposes of illustration, the function 
\code{exp.cov.simple} is provided as
a simple example and implements the R code above. It can also serve as a 
template for creating new covariance functions for the \code{Krig} 
function. }
\seealso{
 Krig, matern.cov, rdist, rdist.earth, gauss.cov, exp.image.cov 
}
\examples{
# exponential covariance matrix ( marginal variance =1) for the ozone
#locations 
out<- exp.cov( ozone$x, theta=100)

# out is a 20X20 matrix
out2<- exp.cov( ozone$x[6:20,],ozone$x[1:2,], theta=100)
# out2 is 15X2 matrix 
# Kriging fit where the nugget variance is found by GCV 
fit<- Krig( ozone$x, ozone$y, exp.cov, theta=100)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
