\name{image.count}
\alias{image.count}
\title{
  Creates image for a 2-d histogram from irregular locations 
}
\description{
Discretizes a matrix of 2-d locations and counts the
number of 
occurrences for each grid box.  
}
\usage{
image.count(x, grid=NULL, nrow=64, ncol=64)
}
\arguments{
\item{x}{

}
\item{grid}{
A list with components x and y of equally spaced values describing the 
centers of the grid points. The default is to use nrow and ncol and the 
ranges of the data locations (x)  to construct a grid. 
}
\item{nrow}{
Number of rows for the image. Note that nrows equals the number of x
values and corresponds to the horizontal dimension.  
}
\item{ncol}{
Number of columns for the image. Note that  
ncol= number of y values and will 
be the vertical dimension.  
}
}
\value{
An list in image format. Components x and y are the grid values , z is a 
nrow X ncol matrix 
with the Z values. NA's are placed at cell locations where Z data has 
not been supplied. 
}
\seealso{
image.smooth, image.plot, Krig.discretize, Krig.replicates  
}
\examples{
look<- image.count( precip$x, nrow=32, ncol=32)
image.plot( look) 
}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
