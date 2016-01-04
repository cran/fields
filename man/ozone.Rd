% fields, Tools for spatial data
% Copyright 2015, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{Chicago ozone test data}
\alias{ChicagoO3}
\alias{ozone}
\title{
  Data set of ozone measurements at 20 Chicago monitoring stations. 
}
\description{
This data set used be named  \code{ozone} but was 
changed to avoid conflict with other packages.
The \code{ChicagoO3} data is a list of 
components, x and y.  x 
component is longitude and latitude position of each of the 20 
Chicago monitoring stations, y is the average 
daily ozone values over the time period 6/3/87-8/30/87.
These data are used extensively for the test scripts and 
simple examples. The lasting scientific value is probably 
minimal. 

}

\format{
This data set is a list containing the following components: 

\describe{
 \item{lon.lat}{ Longitude-latitude positions of monitoring stations.  }
 \item{x}{An approximate Cartesian set of coordinates for the locations 
where the units are in miles. The origin is in the center of the 
locations. }
 \item{y}{ Average daily ozone values over 1987 summer. }
 }
}
\source{
AIRS, the EPA air quality data base. 
}
\seealso{
Tps, Krig  
}
\examples{
fit<- Tps(ChicagoO3$x, ChicagoO3$y) 
# fitting a surface to ozone measurements. 
surface( fit, type="I")
}
\keyword{datasets}
% docclass is function
% Converted by Sd2Rd version 1.21.
