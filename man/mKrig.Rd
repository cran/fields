% fields, Tools for spatial data
% Copyright 2004-2007, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{mKrig}
\alias{mKrig}
\alias{predict.mKrig}
\alias{mKrig.coef}
\alias{mKrig.trace}
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
    m = 2, chol.args=NULL,cov.args=NULL, find.trA = TRUE, NtrA = 20,
                 iseed = 123,  ...)

\method{predict}{mKrig}( object, xnew=NULL,ynew=NULL, derivative=0, ...)

\method{print}{mKrig}( x, digits=4,... )
 mKrig.coef(object, y)

mKrig.trace( object, iseed, NtrA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Matrix of unique spatial locations (or in print or surface 
       the returned mKrig object.) }
  \item{y}{ Vector of observations at spatial locations, missing values
          are not allowed! Or in \code{mKrig.coef} a new vector of 
         observations. }
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
 \item{find.trA}{If TRUE will estimate the effective degrees of freedom using a simple Monte Carlo
method. This will  add to the computational  by approximately \code{NtrA} solutions of the
linear system but the cholesky decomposition is reused.}
\item{NtrA}{Number of Monte Carlo samples for the trace. But if NtrA is greater than or equal to the 
number of observations the trace is computed exactly.}
\item{iseed}{ Random seed ( using \code{set.seed(iseed)}) used to generate iid normals for Monte 
Carlo estimate of the trace.}

  \item{\dots}{In \code{mKrig} and \code{predict} additional arguments that will be passed to the covariance 
                 function. }
  \item{object}{Object returned by mKrig. (Same as "x" in the print function.)}
  \item{xnew}{Locations for predictions.}
   \item{ynew}{New observation vector. \code{mKrig} will reuse matrix decompositions and find the
new fit to these data.}
  \item{derivative}{If zero the surface will be evaluated. If 
    not zero the matrix of partial derivatives will be computed.}
\item{digits}{Number of significant digits used in printed output.}
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

\code{mKrig.coef} finds the "d" and "c" coefficients represent the 
solution using the previous cholesky decomposition for a new data 
vector. This is used in computing the prediction standard error in 
predict.se.mKrig and can also be used to evalute the estimate 
efficiently at new vectors of observations provided the locations and 
covariance remain fixed. 

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

\code{mKrig.trace} This is an internal function called by \code{mKrig}
to estimate the effective degrees of freedom.  The Kriging surface estimate
at the data locations  is
a linear function of the data and can be represented as A(lambda)y.
The trace of A is one useful
measure of the effective degrees of freedom used in the surface
representation. In particular this figures into the GCV estimate of the smoothing parameter. 
It is computationally intensive to find the trace
explicitly but there is a simple Monte Carlo estimate that is often
very useful.  If E is a vector of iid N(0,1) random variables then the
trace of A is the expected value of t(E)AE. Note that AE is simply
predicting a surface at the data location using the synthetic
observation vector E. This is done for \code{NtrA} independent N(0,1)
vectors and the mean and standard deviation are reported in the
\code{mKrig} summary.  Typically as the number of observations is
increased this estimate becomse more accurate. If NtrA is as large as
the number of observations (\code{np}) then the algorithm switches to
finding the trace exactly based on applying A to \code{np} unit
vectors.
}
\value{
  \item{d}{Coefficients of the polynomial fixed part. }
  \item{c}{ Coefficients of the nonparametric part.}
  \item{nt}{ Dimension of fixed part.}
  \item{np}{ Dimension of c.}
  \item{lambda.fixed}{The fixed lambda value}
  \item{x}{Spatial locations used for fitting.}
  \item{knots}{The same as x}
  \item{cov.function.name}{Name of covariance function used.}
 \item{args}{ A list with all the covariance arguments that were specified
 in the call.}
 \item{m}{Order of fixed part polynomial.}
  \item{chol.args}{ A list with all the cholesky arguments that were specified
         in the call.}
 \item{call}{ A copy of the call to mKrig.}
 \item{non.zero.entries}{ Number of nonzero entries in the covariance
matrix for the process at the observation locations.}
 \item{shat.MLE}{MLE of sigma.}
 \item{rho.MLE}{MLE or rho.}
 \item{rhohat}{Estimate for rho adjusted for fixed model degrees of freedom (ala REML).}
 \item{lnProfileLike}{log Profile likelihood for lambda}
  \item{lnDetCov}{Log determinant of the covariance matrix for the observations having factored out
rho.}
\item{Oemga}{GLS covariance for the estimated parameters in the fixed part of the model (d coefficients0.}
\item{qr.VT, Mc}{QR and cholesky matrix decompositions needed to recompute the estimate for new observation vectors.}
\item{fitted.values, residuals}{Usual predictions from fit.}
\item{eff.df}{Estimate of effective degrees of freedom. Either the mean of the Monte Carlo sample or the 
exact value. }
\item{trA.info}{If NtrA ids less than \code{np} then the individual members of the Monte Carlo sample and 
\code{sd(trA.info)/ sqrt(NtrA)} is an estimate of the standard error. If NtrA is greater than or equal to  \code{np} then these are the diagonal elements of A(lamdba).}
\item{GCV}{Estimated value of the GCV function.}
\item{GCV.info}{Monte Carlo sample of GCV functions}
}
\author{Doug Nychka, Reinhard Furrer }
\seealso{Krig, surface.mKrig, Tps, fastTps, \code{ \link{mKrig.grid}}}
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

# interpolate using tapered version of the exponential, 
# the taper scale is set to 1.5 default taper covariance is the Wendland.
# Tapering will done at a scale of 1.5 relative to the scaling 
# done through the theta  passed to the covariance function.


mKrig( x,y,cov.function="stationary.taper.cov",
       theta = 2.0, lambda=.01,
     Taper="Wendland",  Taper.args=list(theta = 1.5, k=2, dimension=2)
           ) -> out2

predict.surface( out2)-> out.p
surface( out.p)

# Try out GCV on a grid of lambda's.
# For this small data set 
# one should really just use Krig or Tps but this is an example of
# approximate GCV that will work for much larger data sets using sparse 
# covariances and the Monte Carlo trace estimate
#
# a grid of lambdas:
lgrid<- 10**seq(-1,1,,15) 

GCV<- matrix( NA, 15,20)
trA<-  matrix( NA, 15,20)
GCV.est<- rep( NA, 15)
eff.df<- rep( NA, 15)
logPL<- rep( NA, 15) 
# loop over lambda's
for (  k in 1:15){
out<- mKrig( x,y,cov.function="stationary.taper.cov",
       theta = 2.0, lambda=lgrid[k],
     Taper="Wendland",  Taper.args=list(theta = 1.5, k=2, dimension=2)  ) 

GCV[k,]<- out$GCV.info
trA[k,]<- out$trA.info
eff.df[k]<- out$eff.df
GCV.est[k]<- out$GCV
logPL[k]<- out$lnProfileLike
}
#
# plot results different curves are for individual estimates and 
# the two lines are whether one averages first the traces or the GCV criterion.
#
par( mar=c(5,4,4,6))
matplot( trA, GCV, type="l", col=1, lty=2, xlab="effective degrees of freedom", ylab="GCV")
lines( eff.df, GCV.est, lwd=2, col=2)
lines( eff.df, rowMeans(GCV), lwd=2)

# add exact GCV computed by Krig 
out0<-  Krig( x,y,cov.function="stationary.taper.cov",
       theta = 2.0, 
     Taper="Wendland",  Taper.args=list(theta = 1.5, k=2, dimension=2), spam.format=FALSE)  
lines( out0$gcv.grid[,2:3], lwd=4, col="darkgreen")

# add profile likelihood 
utemp<- par()$usr
utemp[3:4] <- range( -logPL)
par( usr=utemp)
lines( eff.df, -logPL, lwd=2, col="blue", lty=2)
axis( 4)
mtext( side=4,line=3, "-ln profile likelihood", col="blue")
title( "GCV ( green = exact) and  -ln profile likelihood", cex=2)




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

 hold<- nearest.dist( x, delta= dguess,upper=TRUE)
# exploit spam format to get quick of number of nonzero elements in each row

 hold2<-  diff( hold@rowpointers)
 #  min( hold2) = 55   which we declare close enough 

# now the following will use no less than 55 nearest neighbors 
# due to the tapering. 

\dontrun{
 mKrig( x,y, cov.function="wendland.cov",k=2, theta=dguess, 
            lambda=1e2) ->  look2
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial }
