% fields, Tools for spatial data
% Copyright 2004-2007, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{Colorado Monthly Meteorological Data}
\alias{ CO.elev}
\alias{CO.id}
\alias{CO.loc}
\alias{CO.names}
\alias{CO.ppt}
\alias{CO.ppt.MAM}
\alias{CO.tmax}
\alias{CO.tmax.MAM}
\alias{CO.tmin}
\alias{CO.tmin.MAM}
\alias{CO.years}
\title{Monthly surface meterology for Colorado 1895-1997}

\description{

Source:

 These is a group of R data sets for monthly min/max temperatures and
precipitation over the period 1895-1997. It is a subset extracted from
the more extensive US data record described in at
\url{http://www.image.ucar.edu/Data/US.monthly.met}. Observed monthly
precipitation, min and max temperatures for the conterminous US
1895-1997. See also
 \url{http://www.image.ucar.edu/Data/US.monthly.met/CO.shtml} for an
on line document of this Colorado subset. Temperature is in degrees C and
precipitation is total monthly accumulation in millimeters. Note that
minimum (maximum) monthly tempertuare is the mean of the daily minimum
(maximum) temperatures. 

Data domain:
 
A rectagular lon/lat region [-109.5,-101]x [36.5,41.5] larger than the
boundary of Colorado comprises approximately 400 stations. Although
there are additional stations reported in this domain, stations that
only report preicipitation or only report temperatures have been
excluded. In addition stations that have mismatches between locations
and elevations from the two meta data files have also been excluded. The
net result is 367 stations that have colocated temperatures and
precipitation. 


}

\format{
This group of data sets is organized with the following  objects:
\describe{
\item{CO.info}{A data frame with columns: station id, elev, lon, lat, station name}

\item{CO.elev}{elevation in meters}

\item{CO.id}{ alphanumeric station id codes}

\item{CO.loc}{locations in lon/lat}

\item{CO.ppt CO.tmax CO.tmin}{Monthly means as three dimensional arrays ( Year, Month, Station).
Temperature is in degrees C and precipitation in total monthly
accumulation in millimeters.}

\item{CO.MAM.ppt CO.MAM.tmax CO.MAM.tmin}{Spring seasonal means 
(March, April,May) as two dimensional arrays
   (Year, Station).}
}


}

\examples{

data(COmonthlyMet)

#Spatial plot of 1997 Spring average daily maximum temps
 quilt.plot( CO.loc,CO.tmax.MAM[103,]  )
 US( add=TRUE)
 title( "Recorded MAM max temperatures (1997)")

# min and max temperatures against elevation

matplot( CO.elev, cbind( CO.tmax.MAM[103,], CO.tmin.MAM[103,]),
  pch="o", type="p",
  col=c("red", "blue"), xlab="Elevation (m)", ylab="Temperature (C)")
title("Recorded MAM max (red) and min (blue) temperatures 1997")


} 
\keyword{datasets}
% docclass is data
% Converted by Sd2Rd version 1.21.
