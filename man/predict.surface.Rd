\name{predict.surface}
\alias{predict.surface}
\title{
  Evaluates a fitted function as a surface object  
}
\description{
Evaluates a a fitted model on a 2-D grid keeping any other variables constant.
The resulting object is suitable for use with functions for viewing 3-d
surfaces.  
}
\usage{
predict.surface(out, grid.list=NA, extrap=FALSE, chull.mask, model=NA,
nx=30, ny=30,...)
}
\arguments{
\item{out}{
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
\item{chull.mask}{
Whether to restrict the fitted surface to be on a convex hull, NA's
are assigned to values outside the
convex hull. chull.mask should be a sequence of points defining a convex
hull. Default is to form the convex hull from the observations if this
argument is missing (and extrap is false).  
}
\item{extrap}{
Extrapolation beyond the range of the data. If false function will be
restricted to the convex hull of the observed data or the hull defined
from the points from the argument chull.mask.
 
}
\item{model}{
Model information used to evaluate the fitted surface. ( E.g. a
lambda value for a Tps fit).
 
}
\item{\dots}{
Any other arguments to pass to the predict function associated with the fit
object. 
}
\item{nx}{
Number of grid points in X axis. }
\item{ny}{
Number of grid points in Y axis. }

}
\value{
The usual list components for making contour and perspective plots 
(x,y,z) along  
with the component containing the grid list. 
}
\details{
This function creates the right grid using the grid.list information or
the defaults (make.surface.grid), calls the predict function for the
object with these points
and also adding any extra arguments passed in the ... section, and then 
reforms the results as a surface object ( as.surface). To determine the 
what parts of the prediction grid are in the convex hull of the data the 
function in.poly is used. See the argument inflation in this function is 
used to include a small margin around the outside of the polygon so that 
point on convex hull are included. This funny addition is to prevent 
excluding grid points that fall exactly on the ranges of the data. 
 }
\seealso{
Tps, Krig, predict, grid.list, make.surface.grid, as.surface, surface , 
in.poly
}
\examples{
fit<- Tps( BD[,1:4], BD$lnya)  # fit surface to data 

# evaluate on a grid on two 
# variables holding two fixed  
# default surface and contour plot  
# make grid list
grid.list<- list( KCl="x", MgCl2=mean(BD[,2]),KPO4="y",dNTP=mean(BD[,4]))

out.p<- predict.surface(fit,grid.list, extrap=TRUE) 

surface(out.p) 

set.panel(1,1)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
