\name{predict.se}
\alias{predict.se}
\title{
  Standard errors of predictions 
}
\description{
Calculates the standard error of predictions. This is
usually the fitted object from a function estimate such as from Krig
or Tps. 
}
\usage{
predict.se(object, ...)
}
\arguments{
\item{object}{
A fitted model object of a certain class 
}
\item{\dots}{
Additional arguments to be passed to 
a particular method. e.g. a grid.list or model specification. }
}
\value{
A vector of standard errors for the predicted values. 
}
\details{
This function is generic and will call the appropriate function to calculate 
the standard errors for the object class. The prediction standard error is
for
the estimated function or parameters (a mean value) not for the
prediction of a new observation.
}
\seealso{
predict, predict.surface.se, predict.se.Krig  
}
\examples{
fit<-Tps(ozone$x,ozone$y)
predict.se(fit)                             # std errors of predictions 
xg<- expand.grid(seq(-15,15,,40),seq(-20,20,,40)) 
out<- predict.se(fit,xg)
image.plot( matrix( out, 40))
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
