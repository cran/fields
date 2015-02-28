\name{MLESpatialProcess}
\alias{MLESpatialProcess}
\alias{MLESpatialProcess.fast}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Estimates key covariance parameters for a spatial process.
%%  ~~function to do ... ~~
}
\description{
Maximizes the likelihood to determine the nugget variance (sigma^2), the sill
( rho) and the range (theta) for a spatial process.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MLESpatialProcess(
x, y, theta.grid = NULL, cov.function =
                 "stationary.cov", cov.args = list(Covariance =
                 "Matern", smoothness = 1), ngrid = 10, niter = 15, tol
                 = 0.01, Distance = "rdist", nstep.cv = 50, verbose =
                 FALSE, ...)
                           
MLESpatialProcess.fast(x, y, lambda.start = NULL, theta.start = NULL,
 cov.function = "stationary.cov", cov.args = list(Covariance = "Matern", smoothness = 1),
 relative.tolerance = 0.001, Distance = "rdist", verbose = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix of spatial locations with rows indexing location
    and columns the dimension (e.g. longitude/latitude)
}
  \item{y}{ Spatial observations
}
  \item{theta.grid}{ Grid of theta parameter values to use for grid
  search in maximizing the Likelilood. The defualt is do an initial
  grid search on ngrid points with the range at the 3 an d 97
  quantiles of the pairwise distances.If only two points are passed
  then this is used as the range for a sequence of ngrid points.
}
  \item{cov.function}{The name of the covariance function (See help on Krig for details. )
%%     ~~Describe \code{cov.function} here~~
}
  \item{cov.args}{A list with arguments for the covariance functions. These are usually parameters and other options such as the type of 
  	distance function. 
}
\item{lambda.start}{Starting value for lambda (noise to signal variance ratio, sigma^2/rho).
	If NULL then \code{lambda == .5} is used.}
	
  \item{ngrid}{Number of points in grid search over theta. 

}
\item{nstep.cv}{Number of grid points to use in GCV or REML coarse
	search for optimum.}
\item{theta.start}{Starting value for theta (range parameter).
	If NULL then the median of the pairwise distances is used.}

\item{relative.tolerance}{Convergence criterion used in optim. }


  \item{niter}{
Max number of iterations for the golden section search to maximize over theta.
}
  \item{tol}{
  	Tolerance to declare convergence. 
}
  \item{Distance}{
  	Distance function to use in covariance. 
}

\item{verbose}{If TRUE print out intermediate information for debugging. 
}

  \item{\dots}{ Additional arguments to pass to the Krig function. 
}
}
\details{
	MLESpatialProcess is designed to be a robust but perhaps slow function to 
	maximize the likelihood for a Gaussian spatial process. 
	For fixed range the likelihood is maximzed over the nugget and sill
	parameters using the \code{Krig} function. An outer optimization finds the
	maximum over the range parameter as well. See the help(Krig) for details of the restricted 
	maximum likelihood criterion (REML).
	
	MLESpatialProcess.fast uses the \code{optim}
	function to maximize the  likelihood  computed from the \code{mKrig} function. It is more
	efficient in the computation as it does not find a full eigen decomposition with each new value
	of theta and maximizes the likelihood over theta and lambda simultaneously. 
	
	Note the likelihood can be maximized analytically over the  parameters of the fixed part and with the nugget (sigma) and sill (rho)  reduced to the single parameter lambda= sigma^2/rho. So fixing any other covariance parameters the likelihood is maximzed numerically over lambda and theta. The differences between these two functions is due to the differences between the definition of the
	restricted likelihood used in Krig and the conventional likelihood used in mKrig.
	 

}
\value{
	\code{MLESpatialProcess}:
	A list that includes components:
\code{theta.MLE, rho.MLE, sigma.MLE, lambda.MLE} being the maximum
 likelihood estimates of these
parameters. The component \code{REML.grid} is a two column matrix
 with the
first column being the theta grid and the second column being the
 profiled and restricted  likelihood for that value of theta. Here profile means that
  the likelihood has already been evaluated at the maximum over sigma
   and rho for this value of theta.
	      \code{eval.grid} is a more complete "capture" of the 
	      evaluations being a 
	      6 column matrix with the parameters theta, lambda, sigma,
	       rho, profile likelihood and the effective degrees of
	        freedom. This is just last row of 
	      \code{lambda.est} returned by the core function \code{Krig}
	      
	     \code{MLESpatialProcess.fast} here the returned value is limited because this
	     function isbuilt around calls to \code{mKrig}. Returned value is a list with components:
	     \code{pars}, the MLEs for theta, rho, sigma and lambda,  
		  \code{logLikelihood},values of the log likelihood at the maximum, 
		     \code{eval.grid}, a table with the results from evaluating different combinations of 
		     parameters,

            \code{converge}, convergence flag from optim (0=Successfull) and number of evaluations used to find maximum.
            and \code{call}, the calling arguments. 
            
            }

\author{
Doug Nychka
}

\seealso{
Krig, mKrig.MLE, fastTps.MLE, spatialProcess
}
\examples{
	N<- 100
  set.seed(123)
  x<- matrix(runif(2*N), N,2)
  theta<- .2
  Sigma<-  Matern( rdist(x,x)/theta , smoothness=1.0)
  Sigma.5<- chol( Sigma)
  sigma<- .1
  F.true<- t( Sigma.5) \%*\%  rnorm(N)
  Y<-  F.true +  sigma*rnorm(N)
# find MLE for  sigma rho and theta  smoothness fixed  first
# data set
  obj<- MLESpatialProcess( x,Y)
  obj$pars
  # profile likelihood over theta
  plot(obj$eval.grid[,1], obj$eval.grid[,6], xlab="theta",
  ylab= "log Profile likelihood", type="p" )
  xline( obj$pars["theta"], col="red")
# log likelihood surface over theta and  log lambda
image.plot( obj$logLikelihoodSurface$x,
               obj$logLikelihoodSurface$y, obj$logLikelihoodSurface$z,
               xlab="theta (range)", ylab="log lambda" )
# MLE               
  points( obj$pars[1], log(obj$pars[2]), pch=16, col="magenta", cex=1.2)      

# using "fast" version  
obj.fast<- MLESpatialProcess.fast( x,Y)
obj.fast$pars 
# points where likelihood evaluated:
quilt.plot( log( obj.fast$eval.grid[,1:2] ), obj.fast$eval.grid[,7],
                   xlab="log(theta)",ylab="log(lambda)")

# parameters are slightly different due to the differences of REML and the full likelihood.
                 
# example with a covariate  
\dontrun{
data(COmonthlyMet)
obj2<-  MLESpatialProcess( CO.loc, CO.tmean.MAM.climate)
obj3<-  MLESpatialProcess( CO.loc, CO.tmean.MAM.climate, Z= CO.elev)
ind<- !is.na( CO.tmean.MAM.climate)
obj4<-  MLESpatialProcess.fast( CO.loc[ind,], CO.tmean.MAM.climate[ind],
               Z= CO.elev[ind])
# elevation makes a difference
obj2$pars
obj3$pars
obj4$pars
  }
 \dontrun{
# fits for ozone data
data( ozone2) 	 
NDays<- nrow( ozone2$y)
O3MLE<- matrix( NA, nrow= NDays, ncol=4)
dimnames(O3MLE)<- list(NULL, c("theta",  "lambda", "rho",    "sigma"))
for( day in 1: NDays){
	cat( day, " ")
	O3MLE[day,]<- MLESpatialProcess( ozone2$lon.lat, ozone2$y[day,],
	 Distance="rdist.earth")$pars
}
plot( log(O3MLE[,1]), log(O3MLE[,3]))
}  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
