\name{exp.image.cov}
\alias{exp.image.cov}
\title{
  Exponential,Gaussian and "power" covariance family for 2-d gridded 
locations 
}
\description{
Given two sets of locations on a 2-d grid efficiently multiplies a cross
covariance with a vector. 
}
\usage{
exp.image.cov(ind1, ind2, Y, cov.obj = NULL, setup = F, grid, ...)
}
\arguments{
\item{ind1}{
Matrix of indices for first set of locations this is a two column matrix
where each row is the row/column index of the image element. If missing
the default is to use all grid locations.
}
\item{ind2}{
Matrix of indices for second set of locations. If missing this is taken to
be ind2. If ind1 is missing ind2 is coerced to be all grid locations.
}
\item{Y}{
Vector to multiply by the cross covariance matrix. Y must be the same
locations as those referred to by ind2.
}
\item{cov.obj}{
A list with the information needed to do the multiplication by
convolutions. This is usually found by using the returned list when
setup=T.
}
\item{setup}{
If true do not do the multiplication but just return the covariance object
required by this function. 
}
\item{grid}{
A grid list giving the X and Y grids for the image. (See example below.)
This is only required if setup is true. 
}
\item{\dots}{
Any arguments ( e.g. theta and p) to pass to the exp.cov in setting up the
covariance object. This is only required if setup is true.
}
}
\value{
A vector that is the multiplication of the cross covariance matrix with
the vector Y.
}
\details{
This function was provided to do fast computations for large numbers of
spatial locations and supports the conjugate gradient solution in
krig.image. In doing so the observations can be irregular spaced
but their coordinates must be 2-dimensional and be restricted to grid
points.
(The function as.image will take irregular, continuous coordinates and
overlay a grid on them.) 

Returned value: If ind1 and ind2 are matrices where nrow(ind1)=m and
nrow(ind2)=n then the cross covariance matrix, Sigma is an mXn matrix
(i,j) element is the covariance between the grid locations indexed at
ind1[i,] and ind2[j,]. The returned result is Sigma\%*\%Y. Note that one can
always recover the coordinates themselves by evaluating the grid list at
the indices.  e.g.  cbind( grid\$x[ ind1[,1]], grid\$y[ind1[,2])) will give
the coordinates associated with ind1. Clearly it is better just to work
with ind1!
  
Functional Form: Following the same form as exp.cov for irregular
locations, the covariance is defined as exp( -(D.ij **p)) where D.ij is
the Euclidean distance between x1[i,] and x2[j,] but having first been
scaled by theta. Specifically,

D.ij = sqrt(  sum.k (( x1[i,k] - x2[j,k]) /theta[k])**2 ).

Note that if theta is a scalar then this defines an isotropic covariance
function. 

Implementation: This function does the multiplication on the full
grid efficiently by a 2-d FFT. The irregular pattern in Y is handled by
padding with zeroes and once that multiplication is done only the
appropriate subset is returned. 

As an example assume that the grid is 100X100 let big.Sigma denote the big
covariance matrix among all grid points ( If the parent grid is 100x100
then big.Sigma is 10K by 10K !) Here are the computing steps:

temp<- matrix( 0, 100,100)

temp[ ind2] <- Y

temp2<- big.Sigma\%*\% temp

temp2[ind1]

Notice how much we pad with zeroes or at the end throw away!
Here the matrix multiplication is effected through convolution/FFT tricks
to avoid creating and multiplying big.Sigma explicitly. It is often faster
to multiply the regular grid and throw away the parts we do not need then
to deal directly with the irregular set of locations.

Confusion: In this entire discussion Y is treated as vector. However is
one has complete data then Y can also be interpreted as a image matrix
conformed to correspond to spatial locations. See the last example for
this distinction. 
  
}
\seealso{
  exp.cov, smooth.2d, as.image, krig.image 
}
\examples{
# generate a grid and set of indices based on discretizing the locations
# in the precip dataset
data(precip)
out<-as.image( precip$y, x= precip$x)
ind1<- out$ind
grid<- list( x= out$x, y=out$y)
# there are 750 gridded locations 

# setup to create cov.obj for exponential covariance with range= .25
cov.obj<- exp.image.cov( setup=T, grid=grid, theta=.25) 

result<- exp.image.cov( ind1, ind1, Y= rnorm(750),cov.obj=cov.obj)

# The brute force way would be  
# result<- exp.cov( precip$x, precip$x,theta=.25, C=rnorm(750))

# evaluate the covariance between all grid points and the center of the
#grid

Y<- matrix(0,64,64)
Y[32,32]<- 1
result<- exp.image.cov( Y= Y,cov.obj=cov.obj)
image.plot( matrix( result, 64,64 ))
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
