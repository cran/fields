\name{world}
\alias{world}
\title{
  Plot of the world 
}
\description{
Plots quickly, medium resolution outlines of large land masses and 
bodies of water.  
}
\usage{
world(ylim=c(-90, 90), xlim=c(-180, 180), add=F, ...)
}
\arguments{

\item{ylim}{
range of latitudes 
}
\item{xlim}{
range of longitudes 
}
\item{add}{
If true will add the world map to current plot 
}
\item{\dots}{
These are graphical arguments that are passed to the lines function that 
draws outline.  
}
}
\details{
This function uses the FIELDS dataset world.dat for the coordinates.   
}
\examples{
# Draw map in device color # 2 
world( col=3) 
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
