% fields, Tools for spatial data
% Copyright 2015, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{rdist.earth}
\alias{rdist.earth}
\alias{rdist.earth.vec}

\title{
  Great circle distance matrix or vector
}
\description{
Given two sets of longitude/latitude locations, \code{rdist.earth} computes 
the Great circle (geographic) distance matrix among all pairings and 
\code{rdist.earth.vec} computes a vector of pairwise great circle distances 
between corresponding elements of the input locations using the Haversine 
method and is used in empirical variogram calculations.
}
\usage{
rdist.earth(x1, x2, miles = TRUE, R = NULL)
rdist.earth.vec(x1, x2, miles = TRUE, R = NULL)
}
\arguments{
\item{x1}{
Matrix of first set of lon/lat coordinates first column is the
longitudes
and second is the latitudes. 
}
\item{x2}{
Matrix of second  set of lon/lat coordinates first column is the
longitudes
and second is the latitudes. If missing x1 is used. 
}
\item{miles}{
If true distances are in statute miles if false distances in kilometers. 
}
\item{R}{
Radius to use for sphere to find spherical distances. If NULL the radius
is either in miles or kilometers depending on the values of the miles
argument. If R=1 then distances are of course in radians. 
}
}
\value{
The great circle distance matrix if nrow(x1)=m and nrow(
x2)=n then the returned matrix will be mXn. 
}
\details{
Surprisingly the distance matrix is computed efficiently in R by dot products of the
direction cosines. Thanks to Qing Yang for pointing this out a long time 
ago. }
\author{Doug Nychka, John Paige}
\seealso{
  \link{rdist}, \link{stationary.cov}, \link{fields.rdist.near}
}
\examples{
data(ozone2)
out<- rdist.earth ( ozone2$lon.lat)
#out is a 153X153 distance matrix
upper<-  col(out)> row( out)
# histogram of all pairwise distances. 
hist( out[upper])

#get pairwise distances between first 10 and second 10 lon/lat points
x1 = ozone2$lon.lat[1:10,]
x2 = ozone2$lon.lat[11:20,]
dists = rdist.earth.vec(x1, x2)
print(dists)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
