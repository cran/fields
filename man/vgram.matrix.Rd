\name{vgram.matrix}
\alias{vgram.matrix}
\title{
  Computes a variogram from an image 
}
\description{
Computes a variogram for an image taking into account different directions 
and returning summary information about the differences in each of these 
directions.  
}
\usage{
vgram.matrix(dat, R=5, nsum=1:8, collapse=TRUE)
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
\item{nsum}{
The components returned from applying the FIELDS/FUNFITS stats function. 
Default is to return all 8 statistics, n, mean, s.d. quantiles etc.  
}
\item{collapse}{
If true the default just gives variogram as a function of distance between 
points.  
}
}
\value{
If collapse=T 
  
A list with components:  d, a vector of distances for the differences, 
and vgram, the variogram values. This is the tradiational variogram 
ignoring direction.  
If collapse=F 
A list with components: d, a vector of distances for the differences, 
ind, a two column matrix giving the x and y increment used to compute the 
differences, stats, the result of applying the stats function to each set 
of differences,   and vgram, the mean of the differences 
}
\details{
When collapse=F this function is designed to look for departures from 
isotropy but separating the variogram differences according to 
orientation. For small R this runs efficiently because the differences are 
found by subsetting the image matrix.  
For example the mean associated with the row of ind = ( 2,3)  
will be found from all the differences (X(i,j)- X( i+2,j+3)) 
Here X(.,.) are the values for the spatial field.  In this example  
d= sqr(13) and there will be another entry with the same distance but 
corresponding to the direction (3,2). 
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
out<-vgram.matrix( lennon, collapse=FALSE) #this takes a bit of time 

set.panel(2,1)
plot( out$d, out$vgram) 
plot(out$d, out$vgram.robust) 

# different directions do not appear to have different variances

# matrix plot of orientations.

M<- max( out$ind[,1]) +1 
N<- max( out$ind[,2]) +1
# make up matrix of right size 
temp<- matrix( NA, M,N) 
# fill the matrix 
temp[ out$ind+1]<- out$vgram 
# note values for diagonal elements 
image.plot( (1:M)-1, (1:N)-1, temp, xlab="x difference", 
ylab="y difference") 
set.panel(1,1)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
