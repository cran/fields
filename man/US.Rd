% fields, Tools for spatial data
% Copyright 2004-2013, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{US}
\alias{US}
\title{
  Plot of the US with state boundaries 
}
\description{
Plots quickly, medium resolution outlines of the US with the states and 
bodies of water. A simple wrapper for the map function from the maps package. 
}
\usage{
US( ...)
}
\arguments{
\item{\dots}{
These are the arguments that are passed to the map function
from the maps package.
}
}
\details{
The older version of this function (fields < 6.7.2) used the FIELDS dataset US.dat for the
  coordinates. Currenty this has been switched to use the maps package.    
}
\seealso{
world
}
\examples{
# Draw map in device color # 3 
US( col=3) 
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
