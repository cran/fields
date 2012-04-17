% fields, Tools for spatial data
% Copyright 2004-2011, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{predict.surface}
\alias{predict.surface}
\alias{predict.surface.default}
\alias{predict.surface.derivative}
\alias{predict.surface.derivative.default}
\alias{predict.surface.se}
\alias{predict.surface.se.default}

\title{
  Evaluates a fitted function or its standard errors as a surface object  
}
\description{
Evaluates a a fitted model on a 2-D grid keeping any other variables constant.
The resulting object is suitable for use with functions for viewing 3-d
surfaces.  
}
\usage{
\method{predict}{surface}(object,...)

\method{predict.surface}{default}(object, grid.list = NA, extrap = FALSE, chull.mask =
                 NA, nx = 80, ny = 80, xy = c(1, 2), order.variables="xy",
                  verbose=FALSE,...)

\method{predict}{surface.se}(object,...)

\method{predict.surface.se}{default}(object, grid.list = NA, extrap = FALSE, chull.mask =
                 NA, nx = 80, ny = 80, xy = c(1, 2), 
                 order.variables="xy",verbose=FALSE, ...)

\method{predict}{surface.derivative}(object,...)

\method{predict.surface.derivative}{default}(object, grid.list = NULL, nx = 80, ny = 80, ...) 

}
\arguments{
\item{object}{
An object from fitting a function to data. In FIELDS this is usually a
Krig object. 
}
\item{grid.list}{
A list with as many components as variables describing the surface. 
All components should have a single value except the two that give the 
grid points for evaluation. If the matrix or data frame has column names,  
these must appear in the grid list. See the grid.list help file for more
details. If this is omitted and the fit just depends on two variables the
grid will be made from the ranges of the observed variables. 
 
}
\item{extrap}{
 Extrapolation beyond the range of the data. If \code{FALSE} (the 
default) the predictions will be restricted to the convex hull of the observed 
data or the convex hull defined from the points from the argument chull.mask. 
This function may be slightly faster if this logical is set to 
\code{TRUE} to avoid checking the grid points for membership in the 
convex hull. For more complicated masking a low level creation of a bounding 
polygon and testing for membership with \code{in.poly} may be useful. 

}

\item{chull.mask}{
Whether to restrict the fitted surface to be on a convex hull, NA's
are assigned to values outside the
convex hull. chull.mask should be a sequence of points defining a convex
hull. Default is to form the convex hull from the observations if this
argument is missing (and extrap is false).  
}


\item{nx}{
Number of grid points in X axis. }
\item{ny}{
Number of grid points in Y axis. }

\item{xy}{ A two element vector giving the positions for the "X" and "Y"
variables for the surface. The positions refer to the columns of the x 
matrix used to define the multidimensional surface. This argument is 
provided in lieu of generating the grid list. If a 4 dimensional surface
is fit to data then \code{ xy= c(2,4)} will evaluate a surface using the 
second and fourth variables with  variables 1 and 3 fixed at their median 
values. NOTE: this argument is ignored if a grid.list arguments is 
passed. }

\item{order.variables}{
If "xy" the variables in grid.list are taken in order as "x" then "y". 
If "yx" the roles are reversed.  Suppose a grid.list had components
lat, lon, elevation and one wanted  a lon/lat surface at a fixed 
elevation. Then one would set to "yx" to make "x" lon and "y" lat. 
}

 \item{verbose}{If TRUE prints out some imtermediate results for debugging.}

 \item{\dots}{
Any other arguments to pass to the predict function associated with the fit
object. 
}

}
\value{
The usual list components for making contour and perspective plots 
(x,y,z) along  with labels for the x and y variables. For \code{predict.surface.derivative} the component \code{z} is a three dimensional array  with \code{nx}, \code{ny}, 2. }
\details{
 This function creates the right grid using the grid.list information or
the attribute in xg, calls the predict 
function for the
object with these points and also adding any extra arguments passed in the
... section, and then reforms the results as a surface object
(as.surface). To determine the what parts of the prediction grid are in
the convex hull of the data the function \code{in.poly} is used. The
argument inflation in this function is used to include a small margin
around the outside of the polygon so that point on convex hull are
included. This potentially confusing modification is to prevent excluding
grid points that fall exactly on the ranges of the data. Also note that as written there is no computational savings for evaluting only the convex subset compared to the full grid. 

The derivative function has much fewer arguments and is written to have more clarity in the coding. For example, it assumes a 2-d surface with "x" and "y" in the order that they appear in the fitted object. To mask out the convex hull of the derivatives one could use \code{predict.surface} to get the convex hull and use this for the masking of the deriviative.

 }
\seealso{
Tps, Krig, predict, grid.list, make.surface.grid, as.surface, surface, 
in.poly
}
\examples{
fit<- Tps( BD[,1:4], BD$lnya)  # fit surface to data 

# evaluate fitted surface for  first two 
# variables holding other two fixed at median values

out.p<- predict.surface(fit)
surface(out.p, type="C") 

#
# plot surface for second and fourth variables 
# on specific grid. 

glist<- list( KCL=29.77, MgCl2= seq(3,7,,25), KPO4=32.13, 
                     dNTP=seq( 250,1500,,25))

out.p<- predict.surface(fit, glist)
surface(out.p, type="C")

out.p<- predict.surface.se(fit, glist)
surface(out.p, type="C")

# the unbiquitous ozone data set day 16
  data( ozone2)

  obj<- Tps(ozone2$lon.lat,ozone2$y[16,], m=3)
  lookd<- predict.surface.derivative( obj)
  set.panel( 2,1) 
  image.plot( lookd$z[,,1])
  image.plot( lookd$z[,,2])
 
  
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
