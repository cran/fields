% fields, Tools for spatial data
% Copyright 2015, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{world}
\alias{world}
\alias{world.color}
\alias{in.land.grid}
\alias{world.land}
\title{Plot of the world}
\description{
Plots quickly, medium resolution outlines of large land masses. This
is a simple wrapper for the map function from the maps package. 
}
\usage{
world(...)
world.land(  ...)
world.color( ... )
in.land.grid(...)

}
\arguments{
\item{\dots}{Same arguments used by the \code{map} function from the
maps package.}
}
\details{
See the longstanding \code{maps} package for documentation on this
function.  The functions world.land, world.color and in.land.grid have
been depreciated but can be recovered from versions of fields 6.7.1 or
older.
 }
\seealso{US, in.poly, in.poly.grid}
\examples{

world()
# add the US
US( add=TRUE,col="blue")

world( fill=TRUE) # land filled in black 

## Western Europe
world( xlim=c(-10,18),ylim=c(36,60),fill=TRUE, col="darkgreen",
border="green1")

}
\keyword{hplot}
% docclass is function
