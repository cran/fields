\name{drape.plot}
\alias{drape.color}
\alias{drape.plot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Perspective plot draped with colors in the facets.}
\description{
Function to produce the 
usual wireframe perspective plot with the facets being filled by 
according to a colors. By default the colors are assigned from a
color bar based on the z values. \code{drape.color} is used to create
a color matrix different from the z matrix used for the wireframe.}
\usage{
drape.plot( x,y,z,z2=NULL, col = tim.colors(64), zlim = NULL, 
 add.legend = TRUE, horizontal = TRUE, theta = 30, phi = 20, ...)

drape.color(z, col = tim.colors(64), zlim = NULL,
   transperent.color = "white")

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ grid values for x coordinate (or if x is a list the 
components x y and z are used.)}
  \item{y} { grid values for y coordinate}
  \item{z}{A matrix of z heights }
  \item{z2} { A matrix of z values to use for coloring facets. If NULL 
then z is used for this purpose}
  \item{col}{ A color table for the z values that will be used for draping}
  \item{zlim}{ the z limits default is the range of z. }
  \item{add.legend}{ If true a color strip is added as a legend.}
  \item{horizontal}{ If true color strip is put at bottom of the 
plot, if FALSE it is placed vertically on the right side.  }
  \item{theta}{ x-y rotation angle for perspective. }
  \item{phi}{ z-angle for perspective. }
  \item{transperent.color}{ Color to use when given an NA in z} 
  \item{\dots}{ Other arguments that will be passed to the 
          persp function.  }
}
\value{
\code{drape.plot}
If an assignment is made the projection matrix from persp is returned. 
This information can be used to add additional 3-d features to the plot.
See the \code{persp} help file for an example how to add additional points 
and lines using the \code{trans3d} function and also the example below. 


\code{drape.color} If dim( z) = M,N this function returns an  (M-1)X(N-1) 
matrix where each element is a text string specifying the color. 

}
\details{ 

The legend strip may obscure part of the plot. If so add this 
as another step changing the defaults in image.plot. 

When using drape.color just drop the results into the \code{col} argument
of \code{persp}. Having said written this, no surprises what drape.plot
does it calls \code{drape.color} followed by \code{persp} and finally 
the legend strip is added with \code{image.plot}.  

The color assignments are done by dividing up the zlim interval into
equally spaced bins, adding a very small inflation to the z limits.  The
mean z values, comprising an (M-1)X(N-1) matrix, for each facet are
discretized to the bins. The bin numbers then become the indices used for
the color scale. If zlim is not specified it is the range of the z or z2 
matrix is used to generate the ranges of the color bar. Note that this 
may be different than the range of the mean facets. 


This kind of plot is also supported through the wireframed function in the
lattice package. The advantage of this fields version is that it uses the
standard R graphics functions -- and is written in R code.  

The drape plot is also drawn by the fields \code{surface} function with 
\code{type="P"}. 



}
\author{D. Nychka }

\seealso{ image.plot, persp, plot.surface, surface, lattice, trans3d} 
\examples{

# an obvious choice:
# Mr. R's favorite New  Zealand Volcano!
data( volcano)
M<- nrow( volcano)
N<- ncol( volcano)
x<- seq( 0,1,,M)
y<- seq( 0,1,,N)

drape.plot( x,y,volcano, col=terrain.colors(128))-> pm 

#note the projection has been saved in pm
# add a point marking the summit
max( volcano)-> zsummit
ix<- row( volcano)[volcano==zsummit]
iy<- col( volcano)[volcano==zsummit]
trans3d( x[ix], y[iy],zsummit,pm)-> uv
points( uv, col="magenta", pch="+", cex=2)

# overlay volcano wireframe with gradient in x direction. 

dz<- ( 
     volcano[1:(M-1), 1:(N-1)] - volcano[2:(M), 1:(N-1)] +
     volcano[1:(M-1), 2:(N)] - volcano[2:(M), 2:(N)]  
         )/2

# convert dz to a color scale:
  zlim<- range( c( dz), na.rm=TRUE)
  drape.color( dz, zlim =zlim)-> zcol

# wireframe with these colors 
  persp( volcano, col=zcol, theta=30, phi=20)

# add legend using image.plot function 
  image.plot( zlim=zlim, legend.only =TRUE, horizontal =TRUE)
 

}
\keyword{hplot}% at least one, from doc/KEYWORDS
