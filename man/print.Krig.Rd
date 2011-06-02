% fields, Tools for spatial data
% Copyright 2004-2011, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{print.Krig}
\alias{print.Krig}
\title{
  Print kriging fit results. 
}
\description{
Prints the results from a fitting a spatial process estimate (Krig)
}
\usage{
\method{print}{Krig}(x,digits=4,...)
}
\arguments{
\item{x}{
Object from Krig function. 
}
\item{digits}{
Number of significant digits in printed output. Default is 4. 
}
\item{\dots}{ Other arguments to print.}
}
\value{
Selected summary results from Krig. 
}
\seealso{
print, summary.Krig, Krig 
}
\examples{
fit<- Krig(ozone$x,ozone$y, theta=100) 
print(fit) # print the summary 
fit # this will work too 
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
