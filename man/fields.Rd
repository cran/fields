% fields, Tools for spatial data
% Copyright 2004-2013, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{fields}
\alias{fields-package}
\alias{fields}
\title{
fields - tools for spatial data 
}
\description{
 Fields is a collection of programs for curve and function
 fitting with an emphasis on spatial data and spatial statistics. The
 major methods implemented include cubic and thin plate splines, 
 universal
 Kriging and Kriging for large data sets. One main feature is any
 covariance function implemented in R code can be used for spatial prediction. Another important feature is that fields will take advantage of compactly supported covariance functions in a seamless way through 
the spam package. See \code{library( help=fields)} for a listing of all the
fields contents. 

fields stives to have readable and tutorial code. Take a look at the 
source code for \code{Krig} and \code{mKrig} to see how things work 
"under the hood". 
To load fields with the comments retained in the source
use \code{ keep.source = TRUE} in the \code{library} command. 
We also keep the source on-line:
browse the directory 
\url{http://www.image.ucar.edu/~nychka/Fields/Source} for commented source. 
\url{http://www.image.ucar.edu/~nychka/Fields/Help/00Index.html} is a
page for html formatted help files. (If you obtain the source version of the
package (file ends in .gz) the commented source code is the R subdirectory.) 

\strong{Major methods} 
\itemize{ \item \code{Tps} Thin Plate spline
regression including GCV and REML estimates for the smoothing parameter.  
\item \code{Krig} Spatial process estimation
(Kriging) including support for conditional simulation. 


The Krig function allows you to supply a covariance function that is
written in native R code.  See (\code{stationary.cov}) that includes
several families of covariances and distance metrics including the
Matern and great circle distance. 
\item \code{mKrig} (micro Krig) are
      \code{fastTps}
 fast efficient  Universal Kriging and spline-like functions, that can take advantage of sparse covariance
functions and thus handle very large numbers of spatial locations. 
\item \code{mKrig.MLE} for maximum likelihood estimates of covariance parameters. This function also 
handles replicate fields assumed to be independent realizations at the same locations.  

}

\strong{Other noteworthy functions}
\itemize{

\item \code{vgram} and \code{vgram.matrix} find variograms for spatial data (and
with temporal replications.

\item \code{cover.design} Generates space-filling designs where the distance 
function is expresed in R code.

\item \code{as.image}, \code{image.plot}, \code{drape.plot}, \code{quilt.plot}
\code{add.image}, \code{crop.image}, \code{half.image}, \code{average.image}, 
\code{designer.colors}, \code{color.scale}, \code{in.poly} Many
convenient functions for working with image data and rationally (well,
maybe reasonably) creating and  placing a color scale on an image plot. 
See also  \code{help(grid.list)} for how fields works with grids and \code{US}
and \code{world} for adding a map quickly. 

\item \code{sreg} \code{splint}   Fast 1-D smoothing 
splines and interpolating cubic splines.

}


\strong{ Generic functions that support the methods} 

\code{plot} - diagnostic plots of fit \cr
\code{summary}- statistical summary of fit \cr
\code{print}- shorter version of summary \cr
\code{surface}- graphical display of fitted surface \cr
\code{predict}- evaluation fit at arbitrary points \cr
\code{predictSE}- prediction standard errors at arbitrary points. \cr
\code{sim.rf}- Simulate a random fields on a 2-d grid.

\strong{Getting Started}

 Try some of the examples from help files for \code{Tps} or 
\code{Krig}. 

\strong{Graphics tips}

\code{help( fields.hints)}
 gives some R code tricks for setting up common legends and axes. 
And has little to do with this package!

\strong{Testing}
See \code{help(fields.tests)} for testing fields. 

\strong{Some fields datasets}
\itemize{
\item \code{CO2} Global satelite CO2 concentrations (simulated field)
\item \code{RCMexample} Regional climate model output
\item \code{lennon} Image of John Lennon
\item \code{COmonthlyMet} Monthly mean temperatures and precip for Colorado
\item \code{RMelevation} Digital elevations for the Rocky Mountain Empire
\item \code{ozone2}  Daily max 8 hour ozone concentrations for the US midwest 
for summer 1987.
\item \code{PRISMelevation} Digital elevations for the 
   continental US at approximately 4km resolution
\item \code{NorthAmericanRainfall} 50+ year average and trend for summer rainfall at 
1700+ stations. 
\item \code{rat.diet} Small paired study on rat food intake over time.
\item \code{WorldBankCO2} Demographic and carbon emission data 
 for 75 countries and for 1999. 
}

\strong{DISCLAIMER:}
 The
authors can not guarantee the correctness of any function or program in
this package. 
}
\examples{

# some air quality data, daily surface ozone measurements for the Midwest:
data(ozone2)
x<-ozone2$lon.lat
y<- ozone2$y[16,] # June 18, 1987

# pixel plot of spatial data
quilt.plot( x,y)
US( add=TRUE) # add US map

fit<- Tps(x,y)
# fits a GCV thin plate smoothing spline surface to ozone measurements.
# Hey, it does not get any easier than this!

summary(fit) #diagnostic summary of the fit 

set.panel(2,2)
plot(fit) # four diagnostic plots of  fit and residuals.


set.panel()
surface(fit) # contour/image plot of the fitted surface
US( add=TRUE, col="magenta", lwd=2) # US map overlaid
title("Daily max 8 hour ozone in PPB,  June 18th, 1987")

}
\keyword{datasets}
