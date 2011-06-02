\name{LatticeKrig}
\alias{LatticeKrig}
\alias{MLE.LatticeKrig}
\alias{lattice.precision}
\alias{predict.LatticeKrig}
\alias{print.LatticeKrig}
\alias{surface.LatticeKrig}
\alias{LatticeKrig.coef}
\alias{LatticeKrig.lnPlike}
\title{
Spatial prediction and inference
using a compactly supported multi-resolution basis and a lattice model for 
coefficients.
}
\description{
A variation of dense rank Kriging that uses the Wendland covariance kernel to 
create a regular set of basis functions on a grid. The coefficients of these
basis functions are modeled as a Markov random field. 
}
\usage{
LatticeKrig(x, y = NULL, weights = rep(1, nrow(x)), Z = NULL, NC,
                 lambda = 1, grid.info = NULL, alpha = 1, beta = -0.2,
                 nlevel = 1, eflag = FALSE, iseed = 123, NtrA = 20,
                 setup = FALSE, verbose = FALSE)

MLE.LatticeKrig(x, y, weights = rep(1, length(y)), NC, grid.info =
                 NULL, alpha = 1, beta, lambda, eflag, iseed = 123,
                 NtrA = 20, verbose = FALSE)

lattice.precision(N1,N2, spam.format = TRUE,  beta)
predict.LatticeKrig( object, xnew=NULL,Z=NULL, drop.Z=FALSE, ...)
surface.LatticeKrig( obj, ...)
print.LatticeKrig( x, digits=4, ...)
LatticeKrig.coef(Mc,wPHI,wT.matrix,wy,lambda)
LatticeKrig.lnPlike(Mc, Q, y, lambda, residuals, weights)
}

\arguments{
  \item{x}{Spatial locations}

  \item{y}{Spatial observation}

  \item{weights}{A vector that is proportional to the reciprocal
  variances of the errors. I.e. errors are assumed ot be uncorrealted
  with varainces sigma^2/weights }

  \item{Z}{Linear covariates to be included in fixed part of the model
  that are distinct from the default low order polynomial in \code{x}}
  
  \item{drop.Z}{If true the fixed part will only be evaluated at the
polynomial part of the fixed model. The contribution from the other
covariates will be omitted.}
  
  \item{NC}{Maximum number of grid points in one dimension. If the
domain is rectangular the smaller dimension will have less that NC
points.  }
  
\item{grid.info}{A list with components xmin, xmax, ymin, ymax and
delta that specify ranges and spacing for first level of grid points
that are the centers of the basis functions. If NULL the ranges of the
coordinate locations are used and delta is determined from \code{NC}.}
  
  \item{alpha}{Relative weight given to each level of basis
  functions.}
  
  \item{beta}{The spatial autoregressive parameter that determines the
correlation strength of the Markov random field}
  
  \item{eflag}{If TRUE last resolution level will be tapered
  exponentials.}
  
  \item{nlevel}{Number of levels for the multiresolution basis.}
  
  \item{lambda}{The ratio of the nugget variance to the parameter
  controlling the marginal variance of the process (called rho in
  fields).} \item{iseed}{Random seed used in the Monte Carlo techinque
  for approxmiating the effective degrees of freedom (trace of the
  smoothing matrix) for the GCV criterion. If NA no seed is set. }
  
  \item{NtrA}{ Number of random samples used Monte Carlo method for
determining effective degrees of freedom }
 
  \item{setup}{If FALSE the computations are just setup and no spatial
estimate is computed.  This makes it easy to reuse the symbolic part
of the sparse cholesky decomposition with differents values for the
matrix entries.  See MLE.LatticeKrig to see how this is used to update
the cholesky decomposition.}
 
  \item{spam.format}{If true matrix returned in sparse format.}

  \item{N1}{Number of grid points in x. Precision matrix does not need
the actual centers of the basis functions just the dimensions of the
grid array (N1 and N2).}

  \item{N2}{Number of grid points in y}
  
  \item{object}{The LatticeKrig object}

  \item{obj}{The LatticeKrig object}

  \item{xnew}{Matrix of locations for predition}
 
  \item{digits}{Number of digits in printed output}

  \item{verbose}{If TRUE print out results as they are computed in
  loop.}

  \item{\dots}{Additional optional arguments.}
  
  \item{Mc}{Cholesky decomposition of regression matrix}
 
  \item{Q}{Precision matrix for coefficients.}
  
  \item{wPHI}{weighted matrix of basis functions}
  
  
  \item{wT.matrix}{weighted matrix of fixed part of estimate}

  \item{wy}{weighted observations}

  \item{residuals}{Residuals from fitting spatial process.}

}
\details{
This is an experimental function that uses a specific covariance with the promise of 
providing spatial analysis for very large data sets.  
The model for the spatial field  is 
 
 f(x) = P(x) +  sum Phi.j(x) c.j 

and the data model is 

 Y.k = f(x.k) + e.k

with e.k uncorrelated normal errors with variance sigma**2.

Where P is a linear polynomial (for now), Phi.j  are fixed basis functions and 
c.j are random coefficients. The basis functions are two dimensional radial basis functions (RBF)  based on the Wendland covariance. For a single level (\code{nlevel ==1}) the RBFs
centered at the grid points and with radial support \code{delta*overlap} where \code{delta} is the spacing between grid points and overlap has the default value of 2.5. For multi-resolution each
subsequent level is based on a grid with delta divided by 2. (See example below.)

The prior on c is a multivariate normal, with a mean of zero and the inverse of the covariance matrix, Q, is also known as the  precision matrix. Q is assumed to be block diagonal corresponding to the 
organization of the basis function into levels of resolution.
If \code{nlevels} are specified, the ith block is a precision matrix based on a first order spatial
 autoregression with \code{beta[i]} being the autoregressive parameter. The specific precision matrix for each block (level), Q.i,  is  implemented by \code{lattice.precision}. Briefly, this matrix has 
the form  Q.i= t(H)\%*\%H and is multiplied by the relative weight \code{alpha[i]}. 
Each row of H, corresponding to a point in the lattice in the interior, is one on the diagonal and
 \code{beta[i]} at each of the four nearest neighbors of the lattice point. 
Points on the edges just have less neighbors. The assume is that for the coefficients associated with the ith block of basis functions (c.i),  t(H)\%*\%c.i  will be uncorrelated. 
Values for \code{beta[i]} that are negative and close to 1/4 give reasonable covariance models. The overall design of this function does not require that the centers actually define a grid but then the construction the matrix H may need to modified in the internal function \code{lattice.precision}. 
The reader is referred to the function \code{MR.cov} for an explicit code that computes the implied ccovariance funtion for the process f.

Note that consistent with other fields functions, the two main parameters in the model, sigma**2 and rho are parameterized as  lambda = sigma**2/rho and rho. 

Given this model is possible to compute conditional expectation of f given Y and also the 
profile likelihood for lambda, alpha and beta. Because both the basis functions and Q are sparse it is possible
to compute the estimates and likelihood for large numbers of spatial locations. This model has
an advantage over covariance tapering or compactly supported covariance functions ( e.g. fastTps)
because the covariance functions can have 
longer range correlations. 

The function  \code{MR.cov} is in the form to be used with \code{mKrig} or \code{Krig} and 
is largely for checking and examining the implied covariance function based on the 
grid and the value of \code{nlevel}, \code{alpha}, \code{beta}. 

\code{LatticeKrig} Find spatial process estimate for fixed nlevel, alpha, beta, NC (or delta)
   and lambda.

\code{MLE.LatticeKrig} Simple grid serach over lambda and beta and computes profile likelihood and 
GCV function. 

\code{LatticeKrig.coef} and \code{LatticeKrig.lnPlike} are two low level functions to find the 
basis fucntion coefficients and the evaluate the likelihood. \code{LatticeKrig.lnPlike} also finds the 
profile MLE of sigma and rho given a fixed value for lambda (and alpha, beta, nlevel).
 See the code in LatticeKrig and  also MLE.LatticeKrig to see how these functions are used.

}
\value{
\code{LatticeKrig:} 
A LatticeKrig class object  with components for evaluating the estimate at arbitrary locations, 
describing the fit and returning the symbolic Cholesky decomposition to allow for fast updating with 
new values of lambda and beta. The component \code{MRinfo} is a list with the information that describes the layout of the multiresolution and the basis function can be reconstructed with this list.
(see \code{MRbasis} as an example.)

\code{MLE.LatticeKrig:}
A matrix where each row are the results of evaluating the fit at specific values for beta and lambda.
Columns are: beta, lambda,  trA (Estimated Effective degrees of freedom in fit), SEtrA (Standard of trA estimate), sigmaMLE (MLE for sigma for fixed lambda and beta), 
rhoMLE ( MLE for rho for fixed lambda and beta), lnProfileLike 
(ln likelihood having maximized over rho  and fixed part with lambda and beta fixed.), GCV (GCV function at lambda and beta).



}
\author{
Doug Nychka
}

\seealso{
mKrig, Krig, fastTps, Wendland
}
\examples{

data( ozone2)
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
# LatticeKrig not set up handle missing observations. 

# estimate using a single levellattice with maximum of 20 grid points for x and y.
 obj<- LatticeKrig( x,y,NC=20, lambda=.01, beta=-.2)

# estimate of effective degrees of freedom in surface. 
 obj$trA.est

# plot fitted surface. 
  surface( obj) # see also predict.surface to just evaluate on a grid
  US( add=TRUE)

# list the grid information
  print( obj$MRinfo$grid.info)

# Search over lambda and beta to find MLE.
  beta<- c(-.2,-.21,-.22,-.23,-.24,-.245,-.248, -.249)
  NC<- 20
  lambda<- 10**seq( -7,4,,20)*NC**2
  out<- MLE.LatticeKrig( x,y,NC=NC,alpha=1,beta=beta, lambda=lambda)
#Interpolate these values with Tps
  temp<-Tps( out[,c(2,4)], out[,8],lambda=0)
  surface( temp, xlab="beta", ylab="trA", type="I")
  points(  out[,c(2,4)], pch="o", col="grey50")
  title("ln Profike Likelihood")

# adding in a covariate:

  data(COmonthlyMet)
  y<- CO.tmin.MAM.climate
  good<- !is.na( y)
  y<-y[good]
  x<- as.matrix(CO.loc[good,])
  Z<- CO.elev[good]
  out<- LatticeKrig( x,y,Z=Z,NC=30, lambda=.1, beta=-.25)
  set.panel(2,1)
# quilt.plot with elevations
  quilt.plot( x, predict(out))
# quilt.plot without elevation linear term included
  quilt.plot( x, predict(out, drop.Z=TRUE))
  set.panel()

# a bigger problem:

 data(CO2)
  obj1<- LatticeKrig( CO2$lon.lat,CO2$y,NC=100, lambda=10, beta=-.22)
# 4600 basis functions 100X46 lattice.
  obj1$center.grid.list

  obj1$trA.est # about 1600 effective degrees of freedom 
#
  glist<- list( x= seq( -180,180,,200),y=seq( -80,80,,100) )
  xg<-  make.surface.grid(glist)
  fhat<- predict.LatticeKrig( obj1, xg)
  fhat <- matrix( fhat,200,100) # convert to image
  image.plot( glist$x, glist$y, fhat)
  world( add=TRUE, col="magenta")

# Here is an illustration of using mKrig with this covariance
# to reproduce the computations of LatticeKrig. The difference is that mKrig can
# not take advantage of any sparsity in the precision matrix 
# and the inverse covariance matrix may not be sparse.
#
  beta<- -.2
  lambda <-  1.5
  obj1<- LatticeKrig( x,y,NC=16, lambda=lambda, beta=beta, NtrA=20,iseed=122)
# in both calls iseed is set the same so we can compare 
# Monte Carlo estimates of effective degrees of freedom
  obj2<- mKrig( x,y, lambda=lambda, m=2, cov.function="MR.cov",
                      cov.args=list( grid.info=obj1$MRinfo$grid.info,beta=beta), NtrA=20, iseed=122)
# Example of how tests are automated in fields
  test.for.zero.flag<- TRUE
  test.for.zero( obj1$fitted.values, obj2$fitted.values,
                  tag="comparing predicted values LatticeKrig and mKrig")

# Examing the implied covariance model for this choice.
# grid to evaluate covariance function 
  NP<- 80
  xgrid<- seq( -1,1,,NP)
  xg<-  make.surface.grid(fields.x.to.grid(x))
  grid.info<- obj1$MRinfo$grid.info 
  set.panel(3,1)
# point to evaluate covariance
  x0<- cbind( -106,39) # point close to center
  hold<- MR.cov( xg, x0,grid.info, beta=-.2)
  hold<- as.surface(xg, hold) # reformat as image
  image.plot(hold)

  x0<- cbind( -102,41) # point in corner
  hold<- as.surface(xg,MR.cov( xg, x0,grid.info, beta=-.2))
  image.plot(hold) 
  


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}

