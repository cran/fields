% fields, Tools for spatial data
% Copyright 2004-2011, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{sim.Krig}
\alias{sim.Krig}
\alias{sim.Krig.approx}
\alias{sim.mKrig.approx}
\title{Conditional simulation of a spatial process}
\description{
Generates exact (or approximate) random draws from the conditional 
distribution of a spatial process given specific observations. This is a 
useful way to characterize the uncertainty in the predicted process from 
data. This is known as conditional simulation in geostatistics or 
generating an ensemble prediction in the geosciences. sim.Krig.grid can 
generate a conditional sample for a large regular grid but is restricted 
to stationary correlation functions. 
 }
\usage{
sim.Krig(object, xp, M = 1, verbose = FALSE)

sim.Krig.approx(object, grid.list = NA, M = 1, nx = 40, ny = 40, verbose = 
FALSE, 
 extrap = FALSE)

sim.mKrig.approx(mKrigObject, predictionGrid = NULL,
                            simulationGridlist=NULL, gridRefinement=5, gridExpansion=1,
    M = 1, nx = 40, ny = 40,  verbose = FALSE )
} 
%- maybe also 'usage' for other objects documented here. 
\arguments{
  \item{extrap}{ If FALSE conditional process is not evaluated outside 
     the convex hull of observations. }
  \item{grid.list}{Grid information for evaluating the conditional 
         surface as a grid.list.}
  \item{mKrigObject}{An mKrig Object}
  \item{M}{Number of draws from conditional distribution.}
  \item{nx}{ Number of grid points in x.}
  \item{ny}{ Number of grid points in y.}
  \item{object}{A Krig object.}
  \item{predictionGrid}{A matrix of locations for the conditional simulation.}
  \item{simulationGridlist}{ A \code{gridlist} describing grid for simulation. If missing 
this is created from the range of the locations, \code{nx}, \code{ny} and the other parameters.}
  \item{gridRefinement}{Amount to increase the number of grid points for the simulation grid.}
  \item{gridExpansion}{Amount to increase the size of teh simulation grid. This is used to 
increase the simulation domain so that the circulant embedding algorithm works.}
  \item{xp}{ Locations where to evaluate the conditional process. }
  \item{verbose}{If true prints out intermediate information. }
}

\details{
These functions generate samples from a conditional multivariate 
distribution, or an approximate one, that describes the uncertainty in the estimated spatial 
process under Gaussian assumptions. An important assumption throughout 
these functions is that all covariance parameters are fixed at their 
estimated or prescribed values from the passed object. 

Given a spatial process  Z(x)= P(x) + h(x) observed at 

Y.k = P(x.k) + h(x.k) + e.k

where P(x) is a low order, fixed polynomial and h(x) a Gaussian spatial 
process.
With Y= Y.1, ..., Y.N,
the goal is to sample the conditional distribution of the process. 
 
[Z(x) | Y ]

For fixed a covariance this is just a multivariate normal sampling problem.
\code{sim.Krig.standard} samples this conditional process at the points 
\code{xp} and is exact for  fixed covariance parameters.
\code{sim.Krig.grid} also assumes fixed covariance parameters and does 
approxiamte sampling on a grid. 

The outline of the algorithm is 

0) Find the spatial prediction at the unobserved locations based on the actual data. Call this
 Z.hat(x). 

1) Generate an unconditional spatial process and from this process
simluate synthetic observations.  At this point the approximation is
introduced where the field at the observation locations is
approximated using interpolation from the nearest grid points.

2) Use the spatial prediction model ( using the true covariance) to estimate the 
spatial process at unobserved locations. 

3) Find the difference between the simulated  process and its prediction based on synthetic 
observations. Call this e(x). 

4)  Z.hat(x) + e(x) is a draw from [Z(x) | Y ].

\code{sim.Krig.standard} follows this algorithm exactly. 

\code{sim.Krig.approx} and \code{sim.mKrig.approx} evaluate the conditional surface on grid and
simulates the values of h(x) off the grid using bilinear interpolation
of the four nearest grid points. Because of this approximation it is
important to choose the grid to be fine relative to the spacing of the
observations. The advantage of this approximation is that one can
consider conditional simulation for large grids -- beyond the size
possible with exact methods. Here the method for simulation is
circulant embedding and so is restricted to stationary
fields. The circulant embedding method is known to fail if the domain is small relative to the
correlation range. The argument \code{gridExpansion} can be used to increase the size of the
domain to make the algorithm work. 
}
\value{
\code{sim.Krig}
 a matrix with rowss indexed by the locations in \code{xp} and columns being the 
 \code{M} independent draws.

\code{sim.Krig.approx} a list with components \code{x}, \code{y} and \code{z}. 
x and y define the grid for the simulated field and z is a three dimensional array
 with dimensions \code{c(nx, ny, M)} where the 
first two dimensions index the field and the last dimension indexes the draws.

\code{sim.mKrig.approx} a list with \code{predictionGrid} being the
locations where the field has been simulated. \code{Ensemble} is a
matrix where rows index the simulated vlaues of the field and columns
are the different draws, \code{call} is the calling sequence. Not that if \code{predictionGrid}
has been omitted in the call or is created beforehand using \code{make.surface.grid} it is 
easy to reformat the results into an image format for ploting using \code{as.surface}.
e.g. if \code{simOut} is the output object then to plot the 3rd draw:

\preformatted{
     imageObject<- as.surface(simOut$PredictionGrid, simOut$Ensemble[,3] )
     image.plot( imageObject)
}

}


\author{Doug Nychka}
\seealso{ sim.rf, Krig}
\examples{
data( ozone2)

set.seed( 399)

# fit to day 16 from Midwest ozone data set.
Krig( ozone2$lon.lat, ozone2$y[16,], Covariance="Matern", 
theta=1.0,smoothness=1.0, na.rm=TRUE)-> out

# NOTE theta =1.0 is not the best choice but 
# allows the sim.rf circulant embedding algorithm to 
# work without increasing the domain.

#six missing data locations
 xp<-  ozone2$lon.lat[ is.na(ozone2$y[16,]),]

# 5 draws from process at xp given the data 
# this is an exact calculation
 sim.Krig( out,xp, M=5)-> sim.out

# Compare: stats(sim.out)[3,] to  Exact: predict.se( out, xp)

# simulations on a grid
# NOTE this is approximate due to the bilinear interpolation
# for simulating the unconditional random field. 
# also more  grids points ( nx and  ny) should be used  

sim.Krig.approx(out,M=5, nx=20,ny=20)-> sim.out

# take a look at the ensemble members. 

predict.surface( out, grid= list( x=sim.out$x, y=sim.out$y))-> look

zr<- c( 40, 200)

set.panel( 3,2)
image.plot( look, zlim=zr)
title("mean surface")
for ( k in 1:5){
image( sim.out$x, sim.out$y, sim.out$z[,,k], col=tim.colors(), zlim =zr)
}
\dontrun{
data( ozone2)
y<- ozone2$y[16,]
good<- !is.na( y)
y<-y[good]
x<- ozone2$lon.lat[good,]
O3.fit<- mKrig( x,y, Covariance="Matern", theta=1.0,smoothness=1.0, lambda= .01 )
set.seed(122)
O3.sim<- sim.mKrig.approx( O3.fit, nx=100, ny=100, gridRefinement=4, M=5 )
set.panel(3,2)
surface( O3.fit)
for ( k in 1:5){
image.plot( as.surface( O3.sim$predictionGrid, O3.sim$Ensemble[,k]) )
}
 
}
}
\keyword{spatial}
% at least one, from doc/KEYWORDS
