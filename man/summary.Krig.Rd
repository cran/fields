\name{summary.krig}
\alias{summary.krig}
\title{
  Summary for Krig spatial process estimate 
}
\description{
Creates a list of summary results including estimates for  the nugget
variance (sigma) and the smoothing
parameter (lambda). This list is usually printed using print.summary.Krig.
}
\usage{
summary.Krig(x, digits=4)
}
\arguments{
\item{x}{
A Krig object. 
}
\item{digits}{
Number of significant digits in summary. 
}
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
fit<- Krig(ozone$x, ozone$y, exp.cov, theta=100)  
summary(fit)                            # summary of fit 
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
