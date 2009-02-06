% fields, Tools for spatial data
% Copyright 2004-2007, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{mKrig}
\alias{mKrig}
\alias{predict.mKrig}
\alias{print.mKrig}

\title{"micro Krig"  Spatial process estimate of a curve or surface, 
"kriging" with a known covariance function. }
\description{
This is a simple version of the Krig function that is 
optimized for large data sets and a clear exposition of the
computations. Lambda, the smoothing parameter must be fixed. 
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
  \item{x}{ Matrix of unique spatial locations ( or in print or surface 
       the returned mKrig object.) }
  \item{y}{ Vector of observations at spatial locations, missing values
          are not allowed! }
  \item{weights}{Precision  ( 1/variance) of each observation }
  \item{lambda}{ Smoothing parameter or equivalently the 
ratio between the nugget and process varainces.}
  \item{cov.function}{The name, a text string of the covariance function.}
  \item{m}{The degree of the polynomial used in teh fixed part is 
(m-1) }
  \item{chol.args}{A list of optional arguments (pivot, nnzR)
 that will be used with
the call to the cholesky decomposition. Pivoting is done by default to 
make use of sparse matrices when they are generated. 
This argument is useful in some cases for
sparse covariance functions to reset the memory parameter nnzR.
 (See example below.) }

  \item{cov.args}{A list of optional arguments that will be used in
calls to the covariance function.}

  \item{\dots}{In \code{mKrig} and \code{predict} additional arguments that will be passed to the covariance 
                 function. }
  \item{object}{Object returned by mKrig. (Same as "x" in the print function.)}
  \item{xnew}{Locations for predictions.}
  \item{derivative}{If zero the surface will be evaluated. If 
    not zero the matrix of partial derivatives will be computed.}
}
\details{
 This function is an abridged version of Krig that focuses on the
computations in Krig.engine.fixed done for a fixed lambda parameter
 for unique spatial locations and for data without missing values. These 
restriction simply the code for reading. Note that also little checking 
is done and the spatial locations are not transformed before the 
estimation.  

\code{predict.mKrig} will evaluate the derivatives of the estimated
function if derivatives are supported in the covariance function.
For example the wendland.cov function supports derivatives.

\code{print.mKrig} is a simple summary function for the object.

Sparse matrix methods are handled through overloading the 
usual linear algebra functions with sparse versions. But to take 
advantage of some additional options in the sparse methods the list 
argument chol.args is a device for changing some default values. The 
most important of these is \code{nnzR}, the number of nonzero elements 
anticipated in the Cholesky factorization of the postive definite linear 
system used to solve for the basis coefficients. The sparse of this 
system is essentially the same as the covariance matrix evalauted at the 
observed locations.
As an example of resetting \code{nzR} to 450000 one would use the following 
argument for chol.args in mKrig:

\code{ chol.args=list(pivot=TRUE,memory=list(nnzR= 450000))}

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
\author{Doug Nychka, Reinhard Furrer }
\seealso{Krig, surface.mKrig, Tps, fastTps }
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


# here is a series of examples with a bigger problem 
# using a compactly supported covariance directly

set.seed( 334)
N<- 1000
x<- matrix( 2*(runif(2*N)-.5),ncol=2)
y<- sin( 1.8*pi*x[,1])*sin( 2.5*pi*x[,2]) + rnorm( 1000)*.1
  
mKrig( x,y, cov.function="wendland.cov",k=2, theta=.2, 
            lambda=.1)-> look2

# take a look at fitted surface
predict.surface(look2)-> out.p
surface( out.p)

# this works because the number of nonzero elements within distance theta
# are less than the default maximum allocated size of the 
# sparse covariance matrix. 
#  see  spam.options() for the default values 

# The following will give a warning for theta=.9 because 
# allocation for the  covariance matirx storage is too small. 
# Here theta controls the support of the covariance and so 
# indirectly the  number of nonzero elements in the sparse matrix

\dontrun{
 mKrig( x,y, cov.function="wendland.cov",k=2, theta=.9, lambda=.1)-> look2
}

# The warning resets the memory allocation  for the covariance matirx according the 
# values   'spam.options(nearestdistnnz=c(416052,400))'
# this is inefficient becuase the preliminary pass failed. 

# the following call completes the computation in "one pass"
# without a warning and without having to reallocate more memory. 

spam.options(nearestdistnnz=c(416052,400))
mKrig( x,y, cov.function="wendland.cov",k=2, theta=.9, lambda=1e-2)-> look2

# as a check notice that 
#   print( look2)
# report the number of nonzero elements consistent with the specifc allocation
# increase in spam.options


# new data set of 1500 locations
set.seed( 234)
N<- 1500
x<- matrix( 2*(runif(2*N)-.5),ncol=2)
y<- sin( 1.8*pi*x[,1])*sin( 2.5*pi*x[,2]) + rnorm( N)*.01
  
# the following is an example of where the allocation  (for nnzR) 
# for the cholesky factor is too small. A warning is issued and 
# the allocation is increased by 25% in this example
#
\dontrun{ 
 mKrig( x,y, 
            cov.function="wendland.cov",k=2, theta=.1, 
            lambda=1e2  )-> look2
}
# to avoid the warning 
 mKrig( x,y, 
            cov.function="wendland.cov",k=2, theta=.1, 
            lambda=1e2,
            chol.args=list(pivot=TRUE,memory=list(nnzR= 450000))  )-> look2

# success!

##################################################
# finding a good choice for theta 
##################################################
# Suppose the target is a spatial prediction using roughly 50 nearest neighbors
# (tapering covariances is effective for rouhgly 20 or more in the situation of 
#  interpolation) see Furrer, Genton and Nychka (2006).

# take a look at a random set of 100 points to get idea of scale

set.seed(223)
 ind<- sample( 1:N,100)
 hold<- rdist( x[ind,], x)
 dd<- (apply( hold, 1, sort))[65,]
 dguess<- max(dd)
# dguess is now a reasonable guess at finding cutoff distance for
# 50 or so neighbors

# full distance matrix excluding distances greater than dguess
# but omit the diagonal elements -- we know these are zero!

 hold<- nearest.dist( x, delta= dguess,upper=NULL, diag=FALSE)
# exploit spam format to get quick of number of nonzero elements in each row

 hold2<-  diff( hold@rowpointers)
 #  min( hold2) = 55   which we declare close enough 

# now the following will use no less than 55 nearest neighbors 
# due to the tapering. 

\dontrun{
 mKrig( x,y, cov.function="wendland.cov",k=2, theta=dguess, 
            lambda=1e2) ->  look2
}

#
#    Using mKrig for evaluating  a solution on a big grid.
#    (Thanks to Jan Klennin for motivating this example.)

x<- RMprecip$x
y<- RMprecip$y

Tps( x,y)-> obj

# make up an 80X80 grid that has ranges of observations
# use same coordinate names as the x matrix

glist<- fields.x.to.grid(x, nx=80, ny=80) # this is a cute way to get a default grid that covers x

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
