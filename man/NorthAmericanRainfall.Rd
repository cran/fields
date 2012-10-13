\name{NorthAmericanRainfall}
\alias{NorthAmericanRainfall}
\docType{data}
\title{
Observed North American summer precipitation from the
historical climate network.
}
\description{
Average rainfall in tenths of millimeters for the months
of June, July and August for the period 1950-2010. Data is based on 
1720 stations located in North America.
}
%\usage{data(NorthAmericanRainfall)}
\format{
  The format is a list with components:
"longitude" "latitude"  "precip"    "elevation" "precipSE" "trend" "trendSE" "type"
with elevation in meters, longitude as (-180,180), latitude as (-90, 90) and precipitaion in 1/10 mm 
( precip/254  converts to inches of rainfall)

\code{precip} is the intercept for 1980.5 when a straight line least squares regression is fit to 
each station's record. SE is the companion standard error from the least squares fit.
If the station is complete, then \code{precip} and \code{precipSE} will just be the mean and standard deviation adjusted for a linear trend. The estimated trend  \code{trend} and and its standard error \code{trendSE} are also included. 
 Also due to the centering, for complete data the intercept and trend estimate will be uncorrelated. The component \code{type} indicates whether the station has been "adjusted" (see below) or is still in
  "unadjusted" form. 

}

\source{
The monthly data used to construct this summary was generously provided by Xuebin Zhang, however,
the orignal source is freely available as the Global Historical Climate Network Version 2 Precipitation
quality controlled, curated and served by the US National Climatic Data Center (NCDC). 
The adjusted data from this archive has been modified from its raw form to make the record more homogenous. Heterogenities can come from a variety of sources such as a moving the station a short distance or changes in instruments. See \url{http://www.ncdc.noaa.gov/ghcnm} 
}

\examples{
data(NorthAmericanRainfall)
x<- cbind(NorthAmericanRainfall$longitude,  NorthAmericanRainfall$latitude)
y<- NorthAmericanRainfall$precip
quilt.plot( x,y)
world( add=TRUE)

Zstat<- NorthAmericanRainfall$trend / NorthAmericanRainfall$trendSE
quilt.plot( x, Zstat)

}
\keyword{datasets}
