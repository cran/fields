\name{Wtransform.sim}
\alias{Wtransform.sim}
\title{
  Simulates a 2-d random wavelet field  
}
\description{
Generates a Gaussian field using quadratic W transform basis functions.   
}
\usage{
Wtransform.sim(D, cut.min=16)
}
\arguments{
\item{D}{
A matrix the same size as the image with the variances of the  
wavelet coefficients.   
}
\item{cut.min}{
Coarsest level of resolution in wavelet expansion this set the number of 
smooth basis functions.  
}
}
\value{
A  matrix image of the random field 
}
\section{Side Effects}{

}
\details{
The simulation is easy just an inverse transform applied to weighted 
independent normals. 
The tricky part is getting all the values of D in the right places.  
See plot.coef to plot out the elements of D to check them and also  
Wtransform.D to fill D from variances that are fixed at each level of 
resolution.  
}
\seealso{
Wtransfrom.image, W.image.cov  
}
\examples{
# 
#Fill to look like Gaussian. 
# 
 wght <- c(1., 0.05, 1e-07, 1e-09, 1e-11, 1e-14) 
D<- Wtransform.D(32,32, wght, cut.min=4)$D 
set.panel( 2,2)
for ( k in 1:4){
look<- Wtransform.sim( D, cut.min=4) 
image( look)
}
set.panel( 1,1)

}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
