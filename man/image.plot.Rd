\name{image.plot}
\alias{image.plot}
\title{
  Draws image plot with a legend strip for the color scale.  
}
\description{
This function combines the SPLUS image function with some automatic 
placement of a legend. This is done by splitting the plotting region  
into two parts. Putting the image in one and the legend in the other.  
}
\usage{
image.plot(..., add=TRUE, nlevel=32, legend.shrink=0.9, legend.width=0.02,  
graphics.reset=TRUE, horizontal=FALSE, offset=0.04,  
bigplot=NULL, smallplot=NULL, legend.only=FALSE, col=topo.colors(64))
}
\arguments{
\item{\dots}{
The usual arguments to the image function. 
because this function may change the size of the plotting 
region. (See details below)  
}
\item{add}{
If true add image and a legend strip to the existing plot. 
}
\item{nlevel}{
Number of color levels used in legend strip  
}
\item{legend.shrink}{
 Amount to shrink the size of legend relative to the full height or width
of the plot.  
}
\item{legend.width}{
Width in plotting 
coordinates (the full size plot is [0,1]X[0,1]) of the legend strip. 
}
\item{offset}{
Amount that the legend strip is set in from the left edge or the bottom 
of the plotting region. Units are with respect to the plotting 
coordinates.
}
\item{graphics.reset}{
 If
false (default)  the plotting region ( plt in par) will not be reset and
one can add more information onto the image plot. (e.g. using functions
such as points or lines.)  If true will reset plot parameters to the
values before entering the function.  
}
\item{horizontal}{
If false (default) legend will be a vertical strip on the right side. If
true the legend strip will be along the bottom.  
}
\item{bigplot}{
Plot coordinates for image plot. If not passed
these will
be determined within the function.  
}
\item{smallplot }{
Plot coordinates for legend. If not passed these will be determined within
the function.  
}
\item{legend.only}{
If true just add the
legend to a the plot in the plot region defined by the coordinates in
smallplot.  
}
\item{col}{
Color table to use for image ( see help file on image for details).
Default is a pleasing range of 64 divisions on a topgraphic scale.
}
}
\section{Side Effects}{
After exiting, the
plotting region may be changed to make it possible to add more features to
the plot. To be explicit, par()\$plt may be changed to reflect a smaller
plotting region that includes a legend subplot.  
}
\details{
It is surprising how hard it is just to automatically add the
legend! All "plotting coordinates" mentioned here are in device
coordinates. The plot region is assumed to be [0,1]X[0,1] and plotting
regions are defined as rectangles within this square. We found these
easier to work with than user coordinates.  There are always problems with
default solutions to placing information on graphs but the choices made
here may be useful for most cases. The most annoying thing is that after
using plot.image and adding information the next plot that is made may
have the slightly smaller plotting region set by the image plotting.

 The strategy is simple, divide the plotting region into two smaller
regions.  The image goes in one and the legend in the other. This way
there is always room for the legend. Some adjustments are made to this
rule by not shrinking the image plot if there is already room for the
legend strip and also sticking the legend strip close to the image plot.
Also, one can specify the plot regions explicitly by bigplot and small
plot if the default choices do not work. There may be problems with small
plotting regions in fitting both of these plot and one may have to change
the default character sizes or margins to make things fit.  

 By keeping the zlim argument the same across images one can generate
the same color scale. (See image help file)  One useful technique for a
panel of images is to just draw the first with image.plot to get a legend
and just use image for subsequent plots. Also keep in mind one can just
add a legend to an existing plot without changing plotting parameters.
Usually a square plot (pty="s") done in a rectangular plot region will
have room for the legend with any adjustments stuck to the right side.
}
\seealso{
image
}
\examples{
 x<- 1:10 
 y<- 1:15 
z<- outer( x,y,"+") 
image.plot(x,y,z) 
# now add some points on diagonal  
points( 5:10, 5:10)
#
#fat (5\% of figure) and short (50\% of figure)  legend strip on the bottom 
image.plot( x,y,z,legend.width=.05, legend.shrink=.5, horizontal=TRUE) 

# add a legend on the bottom but first change margin for some room 
par( mar=c(10,5,5,5))

image( x,y,z)
image.plot( zlim=c(0,25), legend.only=TRUE, horizontal=TRUE)
#
 
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
