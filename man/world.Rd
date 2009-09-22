% fields, Tools for spatial data
% Copyright 2004-2007, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{world}
\alias{world}
\alias{world.color}
\title{Plot of the world}
\description{
Plots quickly, medium resolution outlines of large land masses and 
bodies of water.  
}
\usage{
world(ylim = c(-90, 90), xlim = NULL, add = FALSE, asp = 1,
                 xlab = "", ylab = "", xaxt = "n", yaxt = "n", eps =
                 0.1, col=1,shift = FALSE, 
                 fill=FALSE, col.water="white", col.land="darkgrey", ...)

world.color(obj,xlim= c(-180,180), ylim=c(-90,90), 
             col.water="white",col.land="darkgrey", ... )

}
\arguments{
  \item{ylim}{range of latitudes}
  \item{xlim}{range of longitudes}
  \item{add}{logical; if true will add the world map to current plot.}
  \item{asp}{aspect ratio used if \code{add} is false, see
    \code{plot.default} .}
  \item{xlab,ylab}{labels for x- and y-axis; empty by default.}
  \item{xaxt,yaxt}{axis type for x- and y-axis; empty by default.}
\item{eps}{ Tolerance to decide when to insert line break about 0 if
map is to be shifted. (leave this at .1) }
\item{shift}{ If TRUE shifts to be centered on the Dateline and 
longitude runs from 0 to 360. If FALSE centers on Prime Meridian and 
longitude runs 
from -180 to 180.}

\item{col}{Color for map lines when fill is FALSE.}

\item{fill}{ If FALSE draws land outlines. If TRUE fills in land and
water with different colors.}

\item{col.land}{Color for land filling.}

\item{col.water}{Color for water filling.}
\item{obj}{The data set of coast lines \code{world.dat}}
\item{\dots}{
    If the land is not filled these are graphical arguments that are passed to the
    \code{lines} (and \code{plot} if \code{add} is false)
    function that draws the outline.} If fill is TRUE then these arguments
    are passed to the polygon function that does the filling. 

}
\details{
Both functions use the FIELDS dataset world.dat for the coordinates. The 
main advantage of this function is that it is fast and easy to modify. 
The shift option to center over the dateline is useful because often 
plots of oceanic and atmospheric information center the map this way. 

The function \code{world.color} can be used separately but is also
called by world with fill being TRUE. When used alone it will just add
the colored landmasses and water to an existing plot. It is easy to
modify just to add the land masses and use the existing back ground
color as water.  Unfortunately world.color will not work when shift is
TRUE.  The current code could be modified if you need this option. 
Thanks to Steve McIntyre for suggesting and testin the fill option.

 }
\seealso{US}
\examples{

world()
# add the US
US( add=TRUE,col="blue")


world( fill=TRUE) # land filled in dark grey 


## Western Europe (*which* big islands are missing?)
## with a coordinate grid:

world(xlim=c(-10,18),ylim=c(36,60), 
    xaxt = "s", yaxt = "s", fill=TRUE, col.land="darkgreen")

box() # add back in the box that was obscured by the ocean fill.

# add back in outline of land.
world( add=TRUE,lwd=1.5, col="green")

grid()


}
\keyword{hplot}
% docclass is function
