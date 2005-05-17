\name{plot.Krig}
\alias{plot.Krig}
\title{
  Diagnostic and summary plots of the Krig object 
}
\description{
Plots a series of four diagnostic plots that summarize the fit from Krig. 
}
\usage{
plot.Krig(x, main=NA, digits=4, which=c(TRUE,TRUE,TRUE,TRUE), 
graphics.reset=TRUE,
...)
}
\arguments{
\item{x}{
A Krig output object. 
}
\item{main}{
Title of the plot. Default is the function call. 
}
\item{digits}{
Number of significant digits for the RMSE label. 
}
\item{graphics.reset}{
Reset to original graphics parameters after plotting. Default is TRUE.
}
\item{which}{
A vector of 4 logical values. 
Controls which of the four graphs to plot. 
}
\item{\dots}{
Optional graphics arguments to pass to each plot. 
}
}
\details{
This function creates four summary plots of the Krig object. The default
is to put these in a 2X2 panel. However if the screen is already divided
in some other fashion the plots will just be added according to that
scheme. This option is useful to compare to compare several different
model fits.

The first is a scatterplot of predicted value against observed.

The second plot is "standardized" residuals against predicted value.  
Here we mean that the residuals are divided by the GCV estimate for sigma
and multiplied by the square root of any weights that have been specified.
In the case of a "correlation model" the residuals are also divided by the
marginal standard deviation from this model.

The third plot are the values of the GCV function against the effective
degrees of freedom. When there are replicate points several versions of
the GCV function may be plotted.  GCV function is with respect to the
standardized data if a correlation model is specified. A vertical line
indicates the minimium found.

The fourth plot is a histogram of the standardized residuals. 
}
\seealso{
Krig, summary.Krig, Tps  
}
\examples{
fit<-Krig(ozone$x, ozone$y,exp.cov, theta=200)  
# fitting a surface to ozone  
# measurements 
plot(fit)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
