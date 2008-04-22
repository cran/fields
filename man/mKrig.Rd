% fields, Tools for spatial data
% Copyright 2004-2007, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{mKrig}
\alias{mKrig}
\alias{predict.mKrig}
\alias{print.mKrig}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Spatial process estimate of a curve or surface, 
"kriging" with a known covariance function  }
\description{
This is a simple version of the Krig function that is 
optimized for large data sets and a clear exposition of the
computations. 
}
\usage{
mKrig(x, y, weights = rep(1, nrow(x)), 
  lambda = 0, cov.function = "stationary.cov", 
    m = 2, chol.args=NULL,cov.args=NULL, ...)

\method{predict}{mKrig}( object, xnew=NULL, derivative=0, ...)
\method{print}{mKrig}( x, ... )

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ Matrix of unique spatial locations ( or in print the returned 
mKrig object.) }
  \item{y}{ Vector of observations at spatial locations, missing values
          are not allowed! }
  \item{weights}{Precision  ( 1/variance) of each observation }
  \item{lambda}{ Smoothing parameter or equivalently the 
ratio between the nugget and process varainces.}
  \item{cov.function}{ The name, a text string of the covariance function.}
  \item{m}{ The degree of the polynomial used in teh fixed part is 
(m-1) }
  \item{chol.args}{ A list of optional arguments that will be used with
the call to the cholesky decomposition. This is useful in some cases for
sparse covariance functions to reset tmpmax. (See example below.) }

  \item{cov.args}{ A list of optional arguments that will be used in
calls to the covariance function.}

  \item{\dots}{ In \code{mKrig} and \code{predict} additional arguments that will be passed to the covariance 
                 function. }
   \item{object}{ The object returned from mKrig}
   \item{xnew} {Locations for predictions.}
   \item{derivative} { If zero the surface will be evaluated. If 
    not zero the matrix of partial derivatives will be computed.}
}
\details{
 This function is an abridged version of Krig that focuses on the
computations in Krig.engine.fixed done for a fixed lambda parameter
 for unique spatial locations and for data without missing values. These
restriction simply the code for reading. Note that also little checking
is done and the spatial locations are not transformed before the
estimation.  Sparse matrix methods are handled through overloading the
usual linear algebra functions with sparse versions. 

The sparse methods rely on two types of temporary workspace. The
max.points or mean.neighbor arguments used in 
\code{fields.rdist.near} to find nearest
neighboring points  and tmpmax connected with the sparse Cholesky
decompostion. The default values for each of these is a small multiple
of the number of rows of the full matrix. If either of these is too
small reset them to something larger. This is typically the total number
of nonzero elements in the data covariance matrix. 

The temp space for finding near neighbors is 
 max.points = nrow(x1)* mean.neighbor. It is better to set mean.neighbor
rather than max.points because this may then appropriate for different
size matrices that key off of the same distribution of points. This will
be the case in creating the mKrig object and then using this result to
predict to a grid. 

To explicitly reset mean.neighbor as an additional argument to sparse
covariance functions add the argument.

e.g. \code{mean.neighbor =100}


To reset tmpmax include a component in the chol.args list.

e.g. \code{chol.args=list( memory=list( tmpmax=5e5) )}

will reset temp space to 500,000. 

See the example below for more details. 


\code{predict.mKrig} will evaluate the derivatives of the estimated
function if derivatives are supported in the covariance function. 
For example the wendland.cov function supports derivatives.

\code{print.mKrig} is a simple summary function for the object. 

}
\value{
  \item{d}{Coefficients of the polynomial fixed part. }
  \item{c}{ Coefficients of the nonparametric part.}
  \item{nt}{ Dimension of fixed part.}
  \item{np}{ Dimension of c.}
  \item{x}{ Spatial locations used for fitting.}
  \item{cov.function.name}{Name of covariance function used.}
 \item{cov.args}{ A list with all the covariance arguments that were specified
 in the call.}
 \item{chol.args}{ A list with all the cholesky arguments that were specified
         in the call.}
 \item{call}{ A copy of the call to mKrig.}
 \item{non.zero.entries}{ Number of nonzero entries in the covariance
matrix for the process at the observation locations.}
}
\author{Doug Nychka }
\seealso{Krig }
\examples{
#
# Midwest ozone data  'day 16' stripped of missings 
data( ozone2)
y<- ozone2$y[16,]
good<- !is.na( y)
y<-y[good]
x<- ozone2$lon.lat[good,]

# nearly interpolate using defaults (Exponential)
mKrig( x,y, theta = 2.0, lambda=.01)-> out
#
# NOTE this should be identical to 
# Krig( x,y, theta=2.0, lambda=.01) 

# interpolate using tapered version the taper scale is set to 1.5
# Default covariance is the Wendland.
# Tapering will done at a scale of 1.5 relative to the scaling 
# done through the theta  passed to the covariance function.

mKrig( x,y,cov.function="stationary.taper.cov",
       theta = 2.0, lambda=.01, Taper.args=list(theta = 1.5, k=2)
           ) -> out2

predict.surface( out2)-> out.p
surface( out.p)


# here is a bigger problem 
# using a compactly supported covariance directly
# 
# Also note increase in the temp space sie for the
# Cholesky decomposition. Default size (1e5) is too small
#

set.seed( 334)
N<- 1000
x<- matrix( 2*(runif(2*N)-.5),ncol=2)
y<- sin( 1.8*pi*x[,1])*sin( 2.5*pi*x[,2]) + rnorm( 1000)*.01
  
mKrig( x,y, cov.function="wendland.cov",k=2, theta=.1, 
            lambda=1e2)-> look2

# The following will fail for theta=.2 because tmpmax and max.points are too
# small. (Here theta controls the support of the covariance and so 
# indirectly the  number of nonzero elements in the sparse matrix

# mKrig( x,y, cov.function="wendland.cov",k=2, theta=.3, lambda=1e2)-> look2

# as a guess on the size of tmpmax this was set to mean.neighbor* nrow(x)

 mKrig( x,y, 
            cov.function="wendland.cov",k=2, theta=.3, 
            lambda=1e2, mean.neighbor=150, 
            chol.args=list( memory=list( tmpmax=150*1000)) 
             )-> look2



predict.surface( look2)-> out.p
surface( out.p)


#
#    Using mKrig for evaluating  a solution on a big grid.
#    (Thanks to Jan Klennin for motivating this example.)

x<- RMprecip$x
y<- RMprecip$y

Tps( x,y)-> obj

# make up an 80X80 grid that has ranges of observations
# use same coordinate names as the x matrix

glist<- fields.x.to.grid(x, nx=80, ny=80) # this is a cute way to get a defautl grid that covers x

# convert grid list to actual x and y values ( try plot( Bigx, pch="."))
    make.surface.grid(glist)-> Bigx 

# include actual x locations along with grid. 
    Bigx<- rbind( x, Bigx)

# evaluate the surface on this set of points (exactly)

    predict(obj, x= Bigx)-> Bigy

# theta sets range for the compact covariance function 
# this will involve  less than 20 nearest neighbors tha have
# nonzero covariance

    theta<- c( 2.5*(glist$lon[2]-glist$lon[1]), 
                 2.5*(glist$lat[2]-glist$lat[1]))

# this is an interplotation of the values using a compact 
# but thin plate spline like covariance. 
    mKrig( Bigx,Bigy, cov.function="wendland.cov",k=4, theta=theta, 
                 lambda=0)->out2 
# the big evaluation this takes about 45 seconds on a Mac G4 latop
    predict.surface( out2, nx=400, ny=400)-> look

# the nice surface
\dontrun{    
    surface( look)
    US( add=TRUE, col="white")
}


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial }
