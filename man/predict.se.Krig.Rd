\name{predict.se.Krig}
\alias{predict.se.Krig}
\title{
  Standard errors of predictions for Krig spatial process estimate 
}
\description{
Finds the standard error ( or covariance) of prediction based on a linear
combination of
the observed data. The linear combination is usually the "Best Linear
Unbiased Estimate" (BLUE) found from the Kriging equations. There are also
provisions to use a different covariance for evaluation than the one used
to define the BLUE.  
}
\usage{
predict.se.Krig(out, x, cov.function, rho, sigma2, weights=NULL, cov=FALSE,  
stationary=TRUE, fixed.mean=TRUE, ...)
}
\arguments{
\item{out}{
A Krig object.  
}
\item{x}{
Points to compute the SE. 
}
\item{cov.function}{
Name of the true covariance function for the random surface if different from the one  
used to generate predictions. If omitted, then the 
function from the Krig object 
}
\item{rho}{
Parameter that multiplies the covariance function. If omitted the
value 'out\$rho', estimated from the data is used. 
}
\item{sigma2}{
Variance of the measurement error. If omitted the value 'out\$sigma2', 
estimated from the data is used. 
}
\item{weights}{

}
\item{stationary}{
If true the covariance function is assumed to be stationary and more 
efficient computations can be used. 
}
\item{cov}{
If true the full covariance matrix for the predicted values is returned.
Make sure this will not be big if this option is used. ( e.g. 50X50 grid
will return a matrix that is 2500X2500!)  If false just the marginal
standard deviations of the predicted values are returned. Default is
false.
}
\item{fixed.mean}{
If true the mean function is assumed fixed and the variance of this 
component is taken to be zero. If False the variance for this component is 
added to the variance of prediction.   
}
\item{\dots}{
These additional arguments are given to cov.function the "true" covariance 
function, not the one associated with the Kriging weights.   
}
}
\value{
A vector of standard errors for the predicted values of the Kriging fit. 
}
\details{
The predictions are represented as a linear combination of the dependent 
variable, Y. Call this LY.  Based on this representation the conditional 
variance is the same as the expected value of (P(x) + Z(X) - LY)**2.
where 
P(x)+Z(x) is the value of the surface at x and LY is the linear 
combination that estimates this point. Finding this expected value is 
straight forward given the unbiasedness of LY for P(x) and the covariance 
for Z and Y. 

In these calculations it is assumed that the covariance parameters are fixed.
This is an approximation since in most cases they have been estimated from the 
data. It should also be noted that if one assumes a Gaussian field and known 
parameters in the covariance, the usual Kriging estimate is the 
conditional mean of the field given the data. This function finds the 
conditional standard deviations (or full covariance matrix) of the 
fields given the data.  

 There are two useful extensions supported by this function.  Adding
the variance to the estimate of the spatial mean if this is a correlation
model. (See help file for Krig)  and calculating the variances under
covariance misspecification.  Note that the linear combination is based
on the covariance function from the Krig object. One can view this first
step as simply defining a spatial estimator. If the covariance used is
correct it is BLUE, otherwise the MSE for the spatial estimate will be
larger than optimal.  The 'cov.function' argument in this function defaults
to the same covariance used to determine the spatial prediction but
it also can be specified separately, in this case it is interpreted as the true
covariance and the prediction variances are evaluated accordingly. 
}
\section{References}{
See Case Studies in Environmental Statistics  
}
\seealso{
Krig, predict.Krig, predict.surface.se  
}
\examples{
# 
# Note: in these examples predict.se will default to predict.se.Krig using 
# a Krig object  

fit<- Krig(ozone$x,ozone$y,exp.cov, theta=10)    # krig fit 
predict.se.Krig(fit)                        # std errors of predictions 

# make a  grid of X's  
xg<-make.surface.grid( list(seq(-27,34,,40),seq(-20,35,,40)))     
out<- predict.se.Krig(fit,xg)     # std errors of predictions 

#at the grid points out is a vector of length 1600 
# reshape the grid points into a 40X40 matrix etc.  
out.p<-as.surface( xg, out) 
image.plot( out.p) 

# this is equivalent to  the single step function  
# (but default is not to extrapolation beyond data
out<- predict.surface.se( fit) 
image.plot( out) 


# Investigate misspecification 
# 
# first call Krig to create the Krig object.  
# 
Krig( ozone$x, ozone$y, cov.function=exp.cov, theta=100)-> fit 

# note how the new cov. parameters are specified just like in Krig  
predict.se(fit,xg)-> look 
predict.se( fit, xg, cov.function=exp.cov, theta=2.0, sigma2=1)-> look2 

set.panel( 2,1)
image.plot( as.surface( xg, look))
points( fit$x)
image.plot( as.surface( xg, look2))
set.panel( 1,1)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
