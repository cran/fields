% fields, Tools for spatial data
% Copyright 2015, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{summary.Krig}
\alias{summary.Krig}
\alias{summary.spatialProcess}
\title{
  Summary for Krig or spatialProcess estimated models. 
}
\description{
Creates a list of summary results including estimates for  the nugget
variance (sigma) and the smoothing
parameter (lambda). This list is usually printed using  a "print.summary"
function for nice formatting.
}
\usage{
\method{summary}{Krig}(object, digits=4,...)
\method{summary}{spatialProcess}(object, digits=4,...)

}
\arguments{
\item{object}{
A Krig or spatialProcess object. 
}
\item{digits}{
Number of significant digits in summary. 
}
\item{\dots}{Other arguments to summary}
}
\value{
Gives a summary of the Krig object. The components 
include the function call, number of observations, effective degrees  
of freedom, residual degrees of freedom, root mean squared error,  
R-squared and adjusted R-squared, log10(lambda), cost, GCV minimum and  
a summary of the residuals.  
}
\details{
This function is a method for the generic function summary for class
Krig. The results are formatted and printed using print.summary.Krig.
}
\seealso{
Krig, summary, print.summary.Krig  
}
\examples{
fit<- Krig(ChicagoO3$x, ChicagoO3$y, theta=100)  
summary(fit)                            # summary of fit 
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
