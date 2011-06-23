% fields, Tools for spatial data
% Copyright 2004-2011, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{predict.derivative}
\alias{predict.derivative}
\alias{predict.derivative.default}
\title{
Predicted derivatives
}
\description{
Calculates the partial derivatives.
}
\usage{
predict.derivative(object, ...)

\method{predict.derivative}{default}( object, ...)
}
\arguments{
\item{object}{
A fitted model object that support a \code{derivative=1} argument in the
call to predict or a \code{Krig} object. 
}
\item{\dots}{
Additional arguments to be passed usually the locations to evaluate the derivatives.}
}
\value{
The target returned object is intended to be a matrix of parital derivatives.
The columns being the partial derivatives and the rows the locations.
}
\details{
This function is generic with the default to just call \code{predict} with 
\code{derivative=1}. There is a special version for \code{Krig} objects because a separate function was designed to evaluate derivatives. 
}
\seealso{
predict, predict.surface.derivative
}
\examples{
  data(ozone2)
  fit<- Tps( ozone2$lon.lat, ozone2$y[16,]) 
  hold<- predict.derivative(fit)                           
  set.panel(2,1)
  quilt.plot( fit$x, hold[,1])
  quilt.plot( fit$x, hold[,2])
# lazy plots -- better to create a grid of points and evaluate on grid. 
  set.panel()


}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
