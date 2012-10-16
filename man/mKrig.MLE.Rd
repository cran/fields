\name{mKrig.MLE}
\alias{mKrig.MLE}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 maximizes likelihood for the process marginal variance (rho) and nugget standard
 deviation (sigma) parameters (e.g. lambda)
 over a list of covariance models
 or a grid of covariance parameter values.
}
\description{
This function is designed to explore the likelihood surface for different 
covariance parameters having maximized over sigma and rho. This is sometimes 
termed a profile likelihood because sigma and rho are evlauted at there MLE estimates
given the other parameters. 
}
\usage{
mKrig.MLE(x, y, ..., par.grid = NULL, lambda = NULL, lambda.profile = TRUE,  verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{x}{
    Matrix of unique spatial locations (or in print or surface 
    the returned mKrig object.) }
 \item{y}{ 
    Vector or matrix of observations at spatial locations, missing values
    are not allowed! Or in \code{mKrig.coef} a new vector of 
    observations. If y is a matrix the columns are assumed to be 
    independent observations vectors generated from the same covariance and
    measurement error model. }
  \item{\dots}{
     Additional arguments that would also be included in a call to \code{mKrig} to
     specify the covariance model and fixed model covariables.}
  \item{lambda}{If \code{lambda.profile=FALSE} the values of lambda to evaluate the likelihood if \code{TRUE} the starting values for the optimization. If lambda is NA then the optimum value from previous search is used as the starting value. If lambda is NA and it is the first value the starting value defaults to 1.0.} 
     
  \item{par.grid}{
     A list or data frame with components  being
     parameters for  different covariance models. A typical component is \code{theta}
     comprising a vector of scale parameters to try. If par.grid is \code{NULL} then the starting
     then the covariance model is fixed at the components
     in \dots.}
  \item{lambda.profile}{If TRUE maximize likelihood over lambda.}
  \item{verbose}{If TRUE print out interesting intermediate results. }
}
\details{
The observational model follows the same as that described in the \code{Krig} function and 
thus the two primary covariance parameters for a stationary model are the nugget standard deviation
(sigma) and the marginal variance of the process (rho). It is useful to reparametrize as rho and\
lambda= sigma^2/rho. The likelihood can be maximized analytically over 
rho and the parameters in the fixed part of the model the estimate of rho can be substituted back into the likelihood to give a expression that is just a function 
of lambda and the remaining covariance parameters. It is this expression that is then maximized numerically 
over lambda when \code{ lambda.profile = TRUE}.

}
\value{
A list with the components
 \item{out}{A matrix giving the results for evaluating the likelihood for each covariance model.}
\item{par.grid}{The par.grid argument used.}
\item{cov.args.MLE}{The list of covariance arguments (except for lambda) that have the largest 
likelihood over the list covariance models. To fit the surface at the largest likelihood among those tried

\code{ do.call( "mKrig", c(obj$mKrig.args, obj$cov.args.MLE,list(lambda=obj$lambda.opt)) )}

where  \code{obj} is the list returned by this function.}

\item{call}{The calling arguments to this function.}
}
\author{
Douglas W. Nychka
}
\seealso{
Krig.MLE
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# some synthetic data
  N<- 100
  set.seed(123)
  x<- matrix(runif(2*N), N,2)
  theta<- .2
  Sigma<-  Matern( rdist(x,x)/theta , smoothness=1.0)
  Sigma.5<- chol( Sigma)
  sigma<- .1
  M<-5 #  Five (5) independent spatial data sets
  F.true<- t( Sigma.5)\%*\% matrix( rnorm(N*M), N,M)
  Y<-  F.true +  sigma* matrix( rnorm(N*M), N,M)
# find MLE for lambda with range and smoothness fixed in Matern for first data set
  obj<- mKrig.MLE( x,Y[,1], Covariance="Matern", theta=.2, smoothness=1.0)
  obj$out # take a look
  fit<- mKrig( x,Y[,1], Covariance="Matern", theta=.2, smoothness=1.0, lambda= obj$lambda.best) 
#
# search over the range parameter and use all 10 replications for combined likelihood
\dontrun{
  par.grid<- list( theta= seq(.15,.25,,5))
# default starting value for lambda is .02 subsequent ones use previous optimum.
  obj<- mKrig.MLE( x,Y, Covariance="Matern",lambda=c(.02,rep(NA,4)), smoothness=1.0, par.grid=par.grid)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
