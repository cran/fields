% fields, Tools for spatial data
% Copyright 2004-2007, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{vgram.matrix}
\alias{vgram.matrix}
\alias{plot.vgram.matrix}
\title{
  Computes a variogram from an image 
}
\description{
Computes a variogram for an image taking into account different directions 
and returning summary information about the differences in each of these 
directions.  
}
\usage{
vgram.matrix(dat, R=5, dx = 1,dy = 1 )

plot.vgram.matrix(x,...)
}
\arguments{
\item{dat}{
A matrix spacing of rows and columns are assumed to have the same 
distance. 
}
\item{R}{
Maximum radius for finding variogram differences assuming that the grid
points are spaced one unit a part. Default is go out to a 
radius of 5.  
}

\item{dx}{
The spacing of grid points on the X axis. This is used to calculate the 
correct distance between grid points. If dx is not equal to dy then the 
collapse argument must be FALSE. 
}
\item{dy}{ The spacing of grid points on the Y axis.
See additional notes for dx.}

\item{x}{ Returned list from vgram.matrix}

\item{\dots}{ Arguments for image.plot}

}
\value{
  
A list with the following components:  d, a vector of distances for 
the differences, 
and vgram, the variogram values. This is the traditional variogram 
ignoring direction.  

d.full, a vector of distances for all possible shifts up distance R, 
 ind, a two column matrix giving the x and y increment used to compute
the shifts,  and vgram.full, the variogram at each of these 
separations. Also computed is vgram.robust, Cressie's version of a robust 
variogram statistic.

Also returned is the component N the number of differences found for each 
separation csae. 

}

\details{
For the "full" case 
the statistics can summarize  departures from 
isotropy by separating the variogram differences according to 
orientation. For small R this runs efficiently because the differences are 
found by sub-setting the image matrix.  

For example, suppose that a row of the ind matrix is 
(2,3). The variogram value associated with this row is the mean of  
the differences (1/2)*(X(i,j)- X( i+2,j+3))**2 for all i and j. 
(Here X(.,.) are the values for the spatial field.)  In this example  
d= sqrt(13) and there will be another entry with the same distance 
but 
corresponding to the direction (3,2). 
plot.vgram.matrix attempts to organize all the different directions into a 
coherent image plot.
}
\seealso{
\code{\link{vgram}} 
}
\examples{
# variogram for Lennon image.
data(lennon)
out<-vgram.matrix( lennon) 

plot( out$d, out$vgram, xlab="separation distance", ylab="variogram") 
# image plot of vgram values by direction.  

# look at different directions 
out<-vgram.matrix( lennon, R=8)  

plot( out$d, out$vgram) 
# add in different orientations 
points( out$d.full, out$vgram.full, col="red")

#image plot of variogram values for different directions. 
set.panel(1,1)
plot.vgram.matrix( out)
# John Lennon appears remarkably isotropic!

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
