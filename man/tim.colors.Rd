\name{tim.colors}
\alias{tim.colors}
\alias{two.colors}
\title{ Some useful color tables for images. }
\description{
Two color scales useful for image plots:
a pleasing rainbow style color table patterned after that used in 
Matlab by Tim Hoar and also a simple colr interpolation between two 
colors passing through white. 
}
\usage{
tim.colors(n = 64)
two.colors(n=256, start="darkgreen", end="red", middle="white")

}
\arguments{
  \item{n}{ Number of color levels. The setting \code{n}=64 is the orignal 
definition.}
\item{start}{Starting color for lowest values in color scale}
\item{end}{ Ending color.}
\item{middle}{Color scale passes through this color at halfway}
}

\details{
Ask Tim about tim.colors! 
\code{two.colors} is really about three different colors. For other 
colors try \code{fields.color.picker} to view possible choices. 
\code{start="darkgreen", end="azure4"} are the options 
used to get a nice color scale for rendering aerial photos of ski trails. 
(See \url{http://www.image.ucar.edu/Data/MJProject}.)

 
}
\value{
A vector of character strings giving the colors in a hexadecimal format. 

}
\seealso{ topo.colors, terrain.colors, image.plot, quilt.plot, grey.scale }
\examples{

tim.colors(10) 
# returns an array of 10 strings in hex format
#e.g. (red, green,  blue) values of   (16,255, 239)
# translates to "#10FFEF" .


image( outer( 1:20,1:20,"+"), col=tim.colors( 75)) # 75 levels

image( outer( 1:20,1:20,"+"), col=two.colors() )


}
\keyword{ aplot}
