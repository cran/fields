\name{predict.surface.se}
\alias{predict.surface.se}
\title{
  Standard errors of predictions 
}
\description{
Evaluates the standard error of predictions on a surface. 
}
\usage{
predict.surface.se(out, grid.list=NA, extrap=F, chull.mask,...)
}
\arguments{
\item{out}{
A fitted model object of a certain class 
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
restricted to the convex hull of the observed data or the convex hull
defined from the points from the argument chull.mask.
 
}
\item{\dots}{
Any additional arguments that will passed to the predict.se function
specific to the fit object. 
 
}
}
\value{
A surface object with components 
A vector of standard errors for the predicted values. 
}
\details{
This function is generic and will call the appropriate function to calculate 
the standard errors for the object class. It operation is simple  a grid
is created based on the grid list or from the fit object. The prediction
standard error are evaluated on the grid using predict.se. Finally the
standard errors are reformed into a surface object suitable for plotting. 
}
\seealso{
predict, predict.surface, predict.se.Krig, plot.surface, as.surface 
}
\examples{
fit<-Tps(ozone$x,ozone$y)                # tps fit 
out<- predict.surface.se( fit)

surface( out)

# or ... 

image.plot( out)
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
