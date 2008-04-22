% fields, Tools for spatial data
% Copyright 2004-2007, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html
\name{image2lz}
\alias{crop.image}
\alias{half.image}
\alias{get.rectangle}
 \title{Some simple functions for subsetting images} 
\description{These function help in subsetting a image or reducing it
size by averaging adjecent cells.}
\usage{

crop.image(obj, loc=NULL,...)
half.image(obj)
get.rectangle()

}

\arguments{
  \item{obj}{A list in image format with the usual x,y  defining the
grid and z a matrix of image values.}

 \item{loc}{A 2 column matrix of locations within the image region
that define the subset. If not specified then the image is plotted and
the rectangle can be specified interactively.}
\item{\dots}{ Graphics arguments passed to image.plot. This is only relevant
when loc is NULLand the locator function is called. }
}
\details{ 
If loc has more than 2 rows then the largest rectangle
containing the locations is used. 

}

\author{Doug Nychka}

\seealso{ drape.plot, image.plot}

\examples{
data(RMelevation)


loc<- rbind( c(-106.5, 40.8),
             c(-103.9, 37.5))

# extract elevations for just CO frontrange.
FR<- crop.image(RMelevation, loc)
#  check:   image.plot( FR)

# average cells  4 to 1 by doing this twice!
 half.image( RMelevation)-> temp
 half.image( temp)-> temp

# extract the this averaged image
 FR2<- crop.image(temp, loc)
 zr<- range( FR$z)
 set.panel( 1,2)
 image.plot( FR, zlim =zr, horizontal=TRUE)
 image( FR2, zlim =zr, col=tim.colors())
   

}

\keyword{ hplot }% at least one, from doc/KEYWORDS
