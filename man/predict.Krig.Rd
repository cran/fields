\name{predict.Krig}
\alias{predict.Krig}
\title{
  Evaluation of Krig spatial process estimate.  
}
\description{
Provides predictions from the Krig spatial process estimate at arbitrary
points, new data (Y) or other values of the smoothing parameter (lambda)
including a GCV estimate. 
}
\usage{
predict.Krig
(out, x = NULL, lambda = NA, df = NA, model = NA,
        eval.correlation.model = TRUE, y = NULL, verbose = FALSE, gcv = FALSE)
}
\arguments{
\item{out}{
Fit object from the Krig or Tps function. 
}
\item{x}{
Matrix of x values on which to evaluate the kriging surface. 
If omitted, the data x values, i.e. out\$x  will be used. 
}
\item{lambda}{
Smoothing parameter. If omitted, out\$lambda will be used. 
(See also df and gcv arguments)
}
\item{df}{
Effective degrees of freedom for the predicted surface. This can be used
in place of lambda ( see the function Krig.df.to.lambda)
}
\item{model}{
Generic argument that may be used to pass a different lambda.  
}
\item{eval.correlation.model}{
If true ( the default) will multiply the predicted function by marginal
sd's
and add the mean function. This usually what one wants. If false will
return predicted surface in the standardized scale. The main use of this
option is a call from Krig to find MLE's of rho and sigma2     
}
\item{y}{
Evaluate the estimate using the new data vector y (in the same order as
the old data). This is equivalent
to recomputing the Krig object with this new data but is more efficient
because many pieces can be reused. Note that the x values are assumed to
be the same. 
}
\item{verbose}{
Print out all kinds of intermediate stuff for debugging 
}
\item{gcv}{
 
Default is false. Find lambda from the y vector using GCV. See the
returned values for a side effect of setting equal to true. 
}
}
\value{
Vector of predicted responses if gcv=F ( the default)
BUT NOTE: If lambda is found by gcv then the returned values is a list
with the first component a vector of predictions and the second the
estimated value of lambda.  
}
\details{
The main goal in this function is to reuse the Krig object to rapidly
evaluate different estimates. Thus there is flexibility in changing the
value of lambda and also the independent data without having to recompute
the matrices associated with the Krig object. The reason this is possible
is that most on the calculations depend on the observed locations not on
lambda or the observed data.  
}
\seealso{
Krig, predict.surface  
}
\examples{
Krig(ozone$x,ozone$y,exp.cov, theta=50) ->fit
predict( fit) # gives predicted values at data points

grid<- make.surface.grid( list( seq( -40,40,,15), seq( -40,40,,15)))

look<- predict(fit,grid) # evaluate on a grid of points

# some useful graphing functions
out.p<- as.surface( grid, look) # reformat into $x $y $z image-type object
contour( out.p)  


# refit with 10 degrees of freedom in surface

look<- predict(fit,grid, df=15)

# re fit with random data and lambda found by GCV
look<- predict( fit, grid, y= rnorm( 20), gcv=TRUE)

# NOTE: look is a list now  look$predicted  predicted values and look$lambda
# the value of lambda

out.p<-as.surface( grid, look$predicted) 
contour( out.p)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
