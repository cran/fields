\name{REML.test}
\Rdversion{1.1}
\alias{REML.test}
\alias{MLE.Matern}
\alias{MaternGLS.test}
\alias{MaternGLSProfile.test}
\alias{MaternQR.test}
\alias{MaternQRProfile.test}
\title{
Maximum Likelihood estimates for some Matern covariance parameters. 
}
\description{
For a fixed smoothness (shape) parameter these functions provide different ways of
estimating and testing restricted  and profile likehiloods for the Martern covariance
parameters. \code{MLE.Matern} is a simple function that finds the restricted maximum likelihood (REML) estimates of the sill, nugget and range parameters (\code{rho, sigma2 and theta}) of the Matern covariance functions. 
The remaining functions are primarily for testing.
}
\usage{

MLE.Matern(x, y, smoothness, theta.grid = NULL, ngrid=20, verbose=FALSE)

MaternGLSProfile.test(x, y, smoothness = 1.5, init = log(c(0.05,1)))
MaternGLS.test(x, y, smoothness = 1.5, init = log(c(1, 0.2, 0.1)))
MaternQR.test (x, y, smoothness = 1.5, init = log(c(1, 0.2, 0.1))) 
MaternQRProfile.test (x, y, smoothness = 1.5, init = log(c(1))) 

REML.test(x, y, rho, sigma2, theta, nu = 1.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ A matrix of spatial locations with rows indexing location and columns the dimension 
(e.g. longitude/latitude)}
  \item{y}{ Spatial observations
}
  \item{smoothness}{Value of the Matern shape parameter.} 
  \item{theta.grid}{ Grid of theta parameter values to use for  grid search in maximizing the Likelilood. The defualt is do an initial grid search on ngrid points with the range at the 3 an d 97 quantiles of the pairwise distances.If only two points are passed then this is used as the range for a sequence of ngrid points.}
   \item{ngrid}{Number of points in grid search.}
\item{init}{Initial values of the parameters for optimization. For the first three functions these are in the order rho, theta sigma2 and in a log scale. For MaternQRProfile.test initial value is just log(theta). }
\item{verbose}{If TRUE prints more information.}
  \item{rho}{
Marginal variance of Matern process (the "sill")
}
  \item{sigma2}{
Variance of measurement error (the "nugget")
}
  \item{theta}{
Scale parameter (the "range")
}
  \item{nu}{
Smoothness parameter
}

}
\details{ 
 \code{MLE.Matern} is a simple function to find the maximum likelihood 
estimates of using the restricted and profiled likeilihood that is 
intrinsic to the ccomputations in \code{Krig}. The idea is that the 
likelihood is concentrated to the parameters lambda and theta. (where 
lambda = sigma2/rho). For fixed theta then this is maximized over lambda 
using \code{Krig} and thus concetrates the likelihood on theta. The 
final maximization over theta is implemented as a golden section search 
and assumes a convex function. All that is needed is for three theta grid points 
where the middle point has a larger likelihood than the endpoints. In practice the theta grid defualts to a 
20 points equally spaced between the  .03 and .97 quantiles of the distribution of the 
pairwise distances. The likelihood is evaluated at these points and a possible triple is identified. 
If no exists from the grid search the function returns with NAs for the parameter estimates. 
Note that due to the setup of the golden  section search the computation actually minimizes minus the log 
likelihood. 
}
\value{
For MLE.Matern
\item{smoothness}{Value of the smoothness function}
\item{pars}{MLE for rho, theta and sigma}
\item{REML}{Value of  minus the log restricted likelihood at the maxmimum}
\item{trA}{Effective degrees of freedom in the predicted surface based on the MLE parameters.}
\item{REML.grid}{Matrix with values of theta and the log likelihood from the initial grid search.}
}
\author{
Doug Nychka
}
\note{
See the script REMLest.test.R and Likelihood.test.R in the tests directory to 
see how these functions are used to check the likelihood expressions. 
}
\examples{


# Just look at one day from the ozone2 
data(ozone2)

out<- MLE.Matern( ozone2$lon.lat,ozone2$y[16,],1.5)
plot( out$REML.grid)
points( out$pars[2], out$REML, cex=2)
xline( out$pars[2], col="blue", lwd=2)
# to get a finer grid on initial search:
\dontrun{
out<- MLE.Matern( ozone2$lon.lat,ozone2$y[16,],1.5,
                      theta.grid=c(.3,2), ngrid=40) }

# simulated data 
set.seed( 123)
x<- matrix( runif( 200), ncol=2)

rho= 2.0
sigma= .05
theta=.5

Cov.mat<-  rho* Matern( rdist(x,x), smoothness=1.0, range=theta)
A<- chol( Cov.mat)
gtrue<- t(A) \%*\% rnorm( 100)
gtrue<- c( gtrue)
err<-  rnorm(100)*sigma
y<- gtrue + err
out<- MLE.Matern( x,y, 1.0)
# the MLEs:
print( out$pars) 



# Find the MLEs and evaluate the Kriging surface.
data(ozone2)
out<- MLE.Matern( ozone2$lon.lat,ozone2$y[16,],1.5)
#use these parameters to fit surface ....
lambda.MLE<- out$pars[3]/out$pars[1]
out2<- Krig( ozone2$lon.lat,ozone2$y[16,] , Covariance="Matern", theta=out$pars[2],
               smoothness=1.5, lambda= lambda.MLE)
surface( out2) # uses default lambda -- which is the right one.

# here is another way to do this that helps to show how various Krig functions determine the lambda value.
\dontrun{
out2<- Krig( ozone2$lon.lat,ozone2$y[16,] , Covariance="Matern", theta=out$pars[2],
               smoothness=1.5)
# Here the default lambda is that found by GCV
# predict on a grid but use the MLE value for lambda:

out.p<- predict.surface(out2, lambda= lambda.MLE)
surface(out.p) # same surface!
}

# One could also use mKrig with a fixed lambda to compute the surface. 

\dontrun{
# looping  through all the days of the ozone data set.
data( ozone2)
x<- ozone2$lon.lat
y<- ozone2$y
out.pars<- matrix( NA, ncol=3, nrow=89)

for ( k in 1:89){
hold<- MLE.Matern( x,c(y[k,]), 1.5)$pars
cat(k, hold, fill=TRUE)
out.pars[k,]<- hold }
}

}
\keyword{spatial}

