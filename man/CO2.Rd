\name{CO2}
\alias{CO2}
\docType{data}
\title{Simulated global CO2 observations}
\description{
 This is an example of moderately large spatial data set and consist of 
simulated CO2 concentrations.
}
\usage{data(CO2)}
\format{
  The format is a list with two components:
\itemize{
\item lon.lat: 26633x2 matrix of the longitude/latitude locations. 
 These are a subset of  a larger lon/lat grid (see example below).
\item y:   26633 CO2 concentrations in parts per million.
}

}
\details{
This data was provided by Dorit Hammerling and Randy Kawa as a 
test example for the spatial analysis of remotely sensed (i.e. satellite) and 
irregular observations.  
}
\examples{
data(CO2)
#
# A quick look at the observations with world map
quilt.plot( CO2$lon.lat, CO2$y)
world( add=TRUE)

# Note high concentrations in Borneo (biomass burning), Amazonia and
# ... Michigan (???).

# spatial smoothing using the wendland compactly supported covariance
# see help( fastTps) for details
# First smooth using locations and Euclidean distances 
# note taper is in units of degrees 
out<-fastTps( CO2$lon.lat, CO2$y, theta=4, lambda=2.0) 
#summary of fit note about 7300 degrees of freedom 
# associated with fitted surface
 print( out)
# image plot on a grid  (this takes a while)
surface( out, type="I", nx=300, ny=150)
# smooth with respect to great circle distance 
out2<-fastTps( CO2$lon.lat, CO2$y, lon.lat=TRUE,lambda=1.5, theta=4*68) 
print(out2)
#surface( out2, type="I", nx=300, ny=150)

# these data are actually subsampled from a grid. 
# create the image object that holds the data
#
# Note: this code below would not work if some marginal lons or lats
# from the parent grid are missing in the obs. 

x<- unique( CO2$lon.lat[,1])
y<- unique(  CO2$lon.lat[,2])
m<- length( x)
n<- length(y)
z<- matrix( NA, nrow=m, ncol=n)
ind<- cbind( match( CO2$lon.lat[,1], x), match( CO2$lon.lat[,2], y))
z[ind] <- CO2$y

# look at gridded object. 
 image.plot(x,y, z)

# to predict _exactly_ on this grid for the second fit;
# (this take a while)
look<- predict.surface( out2, grid.list=list( x=x, y=y))
image.plot(look)

}
\keyword{datasets}
