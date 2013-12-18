% fields, Tools for spatial data
% Copyright 2004-2013, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{predict.surface}
\alias{predict.surface}
\alias{predict.surface.default}
\alias{predict.surface.derivative}
\alias{predict.surface.se}
\alias{predict.se}
\alias{predict.se.default}
\title{
 This older surface function has been renamed as predictSurface
}
\description{
Evaluates a a fitted model on a 2-D grid keeping any other variables constant.
The resulting object is suitable for use with functions for viewing 3-d
surfaces.  
}
\usage{
\method{predict}{surface}(object,...)

\method{predict.surface}{default}(object, ...)


}
\arguments{
\item{object}{
An object from fitting a function to data. In FIELDS this is usually a
Krig object. 
}

 \item{\dots}{
Any other arguments to pass to the predict function associated with the fit
object. 
}

}
\seealso{
predictSurface, predictSurfaceSE, predictSurfaceDerivative
}

\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
