\name{make.Amatrix}
\alias{make.Amatrix}
\title{
  Matrix relating predicted values to the dependent (Y) values 
}
\description{
For a fixed value of the smoothing parameter or the covariance 
function some nonparametric curve estimates are linear functions of 
the observed data. This is a intermediate  level function that 
computes the linear weights to be applied to the observations to 
estimate the curve at a particular point. For example the predicted 
values can be represented as Ay where A is an N X N matrix of 
coefficients and Y is the vector of observed dependent variables. For 
linear smoothers the matrix A may depend on the smoothing parameter ( 
or covariance function and the independent variables (X) but NOT on Y.    
}
\usage{
make.Amatrix(object, ...)
}
\arguments{
Output object from fitting a data set using a FIELD regression
method. 
Currently this is supported only for Krig ( and Tps) functions. 

\item{object}{

}
\item{\dots}{
Additional arguments that indicate the value of the smoothing 
parameter to use and the X values where predictions should be made.  
}
}
\value{
A matrix where the number of rows is equal to the number of predicted points 
and the number of columns is equal to the length of the Y vector.  
}
\details{
The main use of this function is in finding prediction standard errors. 

For the Krig ( and Tps) functions the A matrix is constructed based on the
representation of the estimate as a generalized ridge regression.  The
matrix expressions are explained in the references from the FIELDS manual.
For linear regression the matrix that gives predicted values is often
referred to as the "hat" matrix and is useful for regression diagnostics.  
For smoothing problems the effective number of parameters in the fit is
usually taken to be the trace of the A matrix. Note that while the A
matrix is usually constructed to predict the estimated curve at the data
point this S function does not have such restrictions. This is possible
because any value of the estimated curve will be a linear function of Y.  
}
\section{References}{
Nychka (2000) "Spatial process estimates as smoothers." 
}
\seealso{
Krig, Tps 
}
\examples{
# Compute the A matrix or "hat" matrix for a thin plate spline 
# check that this gives the same predicted values  
tps.out<-Tps( ozone$x, ozone$y)
A<-make.Amatrix( tps.out, ozone$x)
test<- A\%*\%ozone$y 
# now compare this to predict( tps.out) or tps.out$fitted.values 
#                    they should be the same 
stats( test- tps.out$fitted.values)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
