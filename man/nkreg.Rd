\name{nkreg}
\alias{nkreg}
\title{
  Normal kernel regression estimate 
}
\description{
Estimates a univariate or multivariate regression surface 
using a normal kernel function and a fixed bandwidth. 
}
\usage{
nkreg(data.x, data.y, bandwidth, n.points=50, grid=NULL, grid.list=NULL)
}
\arguments{
\item{data.x }{
a matrix of data locations
}
\item{data.y }{
a vector of values for smoothing 
}
\item{bandwidth}{
bandwidth for kernel if a vector then an estimate will be be found for each 
bandwidth. 
}
\item{n.points }{
Number of equally spaced points to evaluate a univariate estimate 
estimate or a 2-d estimate on a regular grid. 
}
\item{grid}{
matrix of locations to evaluate the kernel estimate
}
\item{grid.list }{
A grid.list that describes the regular grid to evaluate the estimate.
If it is missing the default for 2-d is to create a grid 
}
}
\value{

\item{x}{
Points for evaluation
}
\item{y}{
estimate of curve or surface. If more than one bandwidth is supplied then 
y is matrix with columns indexed by the bandwidth values
}
\item{h}{
bandwidths used. 
}
\item{grid.list}{
grid.list that was either passed or created. 
}
}
\seealso{
Tps, Krig, sreg, smooth.2d, image.smooth
}
\examples{
out<- nkreg( rat.diet$t, rat.diet$con, bandwidth=3, n.points=100)
plot( rat.diet$t, rat.diet$con)
lines( out$x, out$y)

# 2-d example evaluate at 40 points
out<- nkreg( precip$x, precip$y, bandwidth=.5, n.points=64)
image.plot( as.surface(out$x, out$y))
}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
