\name{fields}
\alias{fields}
\title{
fields - tools for spatial data 
}
\description{
Fields is a collection of programs for curve and function
 fitting with an emphasis on spatial data and spatial statistics. The
 major methods implemented include cubic and thin plate splines, universal
 Kriging and Kriging for large data sets. One main feature is any
 covariance function implemented in R can be used for spatial prediction. 

Some major methods include:
\itemize{
\item \code{Tps}   Thin Plate spline regression (including GCV) 
\item \code{Krig}  Spatial process estimation (Kriging) 
\item \code{krig.image}   Spatial process estimate for large problems
}
The Kriging functions allow you to supply a covariance function that is
written in native S/R code. The Matern family (\code{Matern.cov}) is 
included as a built in 
covariance. 

Some other noteworthy functions are
\itemize{
\item \code{cover.design}  Space-filling designs where the distance 
function is expresed in R/S code

\item \code{as.image}    \code{image.plot}  
\code{image.count} \code{smooth.image}  convenient functions 
for working with image data

\item \code{sreg},  \code{qsreg} \code{splint}   Fast 1-D smoothing 
splines and 1-D 
quantile/robust and interpolating cubic splines
}

There are also generic functions that support 
these methods such as 

\code{plot} - diagnostic plots of fit \cr
\code{summary}- statistical summary of fit \cr
\code{print}- shorter version of summary \cr
\code{surface}- graphical display of fitted surface \cr
\code{predict}- evaluation fit at arbitrary points \cr
\code{predict.se}- prediction standard errors at arbitrary points. \cr
\code{sim.rf}- Simulate a random fields on a 2-d grid.

To get started, try some of the examples from help files for \code{Tps} or 
\code{Krig}. See also the manual/tutorial at 
\url{http://www.cgd.ucar.edu/stats/Software/fields}


DISCLAIMER:

This is software for statistical research and not for commercial uses. The
authors do not guarantee the correctness of any function or program in
this package. Any changes to the software should not be made without the
authors permission.
}
\examples{
set.panel()
fit<- Tps(ozone$x, ozone$y)  # fits a thin plate spline surface to ozone 
#                               measurements.
summary(fit) #summary of the fit 
plot(fit) # diagnostic plots of  fit and residuals.
surface( fit, type="I") # image and contour plot of the fitted surface

}
\keyword{datasets}
