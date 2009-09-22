% fields, Tools for spatial data
% Copyright 2004-2007, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{predict.se.Krig}
\alias{predict.se.Krig}
\alias{predict.se.mKrig}
\alias{predict.se.KrigA}
\title{
  Standard errors of predictions for Krig spatial process estimate 
}
\description{
Finds the standard error ( or covariance) of prediction based on a linear
combination of
the observed data. The linear combination is usually the "Best Linear
Unbiased Estimate" (BLUE) found from the Kriging equations. 
This statistical computation is done under the assumption that the
covariance function is known.
}
\usage{
\method{predict.se}{Krig}(object, x = NULL, cov = FALSE, verbose = FALSE,...)
\method{predict.se}{mKrig}(object, x = NULL, verbose = FALSE,...)
 }
\arguments{
\item{object}{
A Krig or mKrig object.  
}
\item{x}{
Points to compute the predict standard error or the prediction 
cross covariance matrix.  
}

\item{cov}{
If TRUE the full covariance matrix for the predicted values is returned.
Make sure this will not be big if this option is used. ( e.g. 50X50 grid
will return a matrix that is 2500X2500!)  If FALSE just the marginal
standard deviations of the predicted values are returned. Default is
FALSE -- of course.
}

\item{verbose}{If TRUE will print out various information for debugging.}

\item{\dots}{
These additional arguments passed to the predict.se function.
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

 There are two useful extensions supported by this function.  Adding the 
variance to the estimate of the spatial mean if this is a correlation 
model. (See help file for Krig) and calculating the variances under 
covariance misspecification.  The function \code{predict.se.KrigA} uses 
the smoother matrix ( A(lambda) ) to find the standard errors or 
covariances directly from the linear combination of the spatial 
predictor. Currently this is also the calculation in 
\code{predict.se.Krig} although a shortcut is used 
\code{predict.se.mKrig} for mKrig objects and this shortcut is planned 
for in a later version of fields

}
\seealso{
Krig, predict.Krig, predict.surface.se  
}
\examples{
# 
# Note: in these examples predict.se will default to predict.se.Krig using 
# a Krig object  

  fit<- Krig(ozone$x,ozone$y,cov.function="Exp.cov", theta=10)    # Krig fit 
  predict.se.Krig(fit)      # std errors of predictions at obs.

# make a  grid of X's  
  xg<-make.surface.grid( 
  list(East.West=seq(-27,34,,20),North.South=seq(-20,35,,20)))     
  out<- predict.se.Krig(fit,xg)   # std errors of predictions 

#at the grid points out is a vector of length 400 
#reshape the grid points into a 20X20 matrix etc.  

   out.p<-as.surface( xg, out) 
   surface( out.p, type="C") 

# this is equivalent to  the single step function  
# (but default is not to extrapolation beyond data
# out<- predict.surface.se( fit) 
# image.plot( out) 



}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
