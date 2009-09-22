% fields, Tools for spatial data
% Copyright 2004-2007, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{RMprecip}
\alias{RMprecip}
\alias{RMelevation}
\alias{PRISMelevation}
\title{
Monthly total precipitation (mm) for August 1963 in the Rocky Mountain 
Region and some gridded 4km elevation data sets.  
}
\description{
\code{RMprecip} is a useful spatial data set of moderate size consisting of 865 
locations. See www.image.ucar.edu/Data for the source of these data. 
\code{PRISMelevation} and \code{RMelevation} are gridded elevations for the 
continental US and Rocky Mountain region at 4km resolution. 
Note that the gridded elevations from the PRISM data product are
different than the exact station elevations. (See example below.)
}

\format{
The data set  \code{RMprecip} is a list containing the following components: 

\describe{
\item{x}{
Longitude-latitude position of monitoring stations. Rows names are station id codes. 
}
\item{elev}{
Station elevation in meters. 
}
\item{y}{
Monthly total precipitation in millimeters.
for August 1963
}
}

The data sets 
 \code{PRISMelevation} and  
 \code{RMelevation} are lists 
in the usual R grid format for images and contouring

They have  the following components:
\describe{
\item{x}{
Longitude-latitude grid at approximately 4km resolution}
\item{elev}{
Average  elevation for grid cell in meters
}

}

These elevations and the companion grid formed the basis for the
103-Year High-Resolution Precipitation Climate Data Set for the
Conterminous United States
\url{ftp://ftp.ncdc.noaa.gov/pub/data/prism100} archived at the National
Climate Data Center. This work was primarily authored by Chris Daly
\url{www.prism.oregonstate.edu} and his PRISM group but had some
contribution from the Geophysical Statistics Project at NCAR. 
and is an interpolation of the observational data to a 4km grid that
takes into account topography such as elevation and aspect.

}

\examples{
# this data set was created  the 
# historical data  taken from 
# Observed monthly precipitation, min and max temperatures for the coterminous US 
# 1895-1997
# NCAR_pinfill 
# see the Geophysical Statistics Project datasets page for the supporting functions 
# and details. 

# plot 
quilt.plot(RMprecip$x, RMprecip$y)
US( add=TRUE, col=2, lty=2)

# comparison of station elevations with PRISM gridded values

data(RMelevation)

interp.surface( RMelevation, RMprecip$x)-> test.elev

plot( RMprecip$elev, test.elev, xlab="Station elevation", 
ylab="Interpolation from PRISM grid")
abline( 0,1,col="blue")

# some differences  with high elevations probably due to complex
# topography!

#
# view of Rockies looking from theSoutheast

save.par<- par(no.readonly=TRUE)

par( mar=c(0,0,0,0))

# fancy use of persp with shading and lighting.
persp( RMelevation, theta=75, phi= 15, 
          box=FALSE, axes=FALSE, xlab="", ylab="", 
         border=NA,
         shade=.95, lphi= 10, ltheta=80,
         col= "wheat4", 
         scale=FALSE, expand=.00025)

# reset graphics parameters and a more conventional image plot.
par( save.par)
image.plot(RMelevation, col=topo.colors(256))
US( add=TRUE, col="grey", lwd=2)
title("PRISM elevations (m)")

} 
\keyword{datasets}
% docclass is data
% Converted by Sd2Rd version 1.21.
