\name{make.Amatrix.Krig}
\alias{make.Amatrix.Krig}
\title{
  Computes the  prediction matrix for a Krig fit. 
}
\description{
For fixed covariance parameters the Krig estimate is a linear function of
the data. This function return the matrix that maps the observed y vector
into the predicted values.
}
\usage{
make.Amatrix.Krig(out, x0=out$x, lambda)
}
\arguments{
\item{out}{
The output object from using Krig to fit a data set 
}
\item{x0}{
The points where the function is to be predicted 
}
\item{lambda}{
Value of the smoothing parameter. Default is value from the Krig
object. 
}
}
\value{
A matrix such that when multiplied times the data vector Y gives the
predicted values 
of the spatial process estimate at the points x0. 
}
\seealso{
Krig, make.Amatrix, predict.se.Krig  
}
\examples{
# Compute the A matrix or "hat" matrix for a spatial process estimate
# with an exponential covariance, range= 100.
# check that this gives the same predicted values 
krig.out<- Krig( ozone$x, ozone$y, exp.cov, theta=100)
A<- make.Amatrix( krig.out, ozone$x)
test.fitted.values<- A\%*\%ozone$y  

# now compare this to predict( krig.out) or krig.out$fitted.values 
#                    they are all  the same! 
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
