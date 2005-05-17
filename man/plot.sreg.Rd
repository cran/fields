\name{plot.sreg}
\alias{plot.sreg}
\title{
  Plots smoothing spline regression object 
}
\description{
Plots a series of three plots that summarizes the smoothing spline regression 
object. 
}
\usage{
plot.sreg(x, main = NA, digits = 4, which = c(TRUE, TRUE, TRUE, TRUE), 
               graphics.reset= TRUE, ...)
}
\arguments{
\item{x}{
A sreg object. 
}
\item{main}{
Title of the plot. Default is the function call. 
}
\item{digits}{
Number of significant digits for the RMSE label. 
}
\item{which}{
Logical vector of length 4 that controls which of the four graphs to
plot. Default is to plot all 4. 
}
\item{graphics.reset}{
Reset to original graphics parameters after plotting. Default is TRUE.
}
\item{\dots}{
Optional plotting arguments that are passed to the first two scatterplots. 
e.g pch="."
}
}
\details{
This function creates four summary plots of the sreg object-- the results 
of fitting a cubic smoothing spline to a data set. The
first plot is a plot of the sreg fit(s) to the data. The second is
predicted values vs 
residuals, and the third is effective number of parameters vs estimated  
average squared prediction errors. The fourth is a histogram ( or
boxplots) of residuals. 
}
\seealso{
sreg, plot  
}
\examples{
fit<-sreg(rat.diet$t,rat.diet$con)
# fit rat data
plot(fit)                       # plot sreg fit 
}
\keyword{smooth}
% docclass is function
% Converted by Sd2Rd version 1.21.
