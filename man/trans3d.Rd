\name{trans3d}
\alias{trans3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Converts 3-d coordinates to 2-d plotting coordinates }
\description{
Uses the projection matrix returned by persp to transform x,y,z values to 
2d plotting region. 
}
\usage{
trans3d(x, y, z, pmat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ vector of x values }
  \item{y}{ vector of y values }
  \item{z}{ vector of z values  }
  \item{pmat}{ projection matrix returned by persp }
}
}
\value{
   A list
  \item{x }{ values of x coordinate in plotting scale}
  \item{comp2 }{values of y coordinate in plotting scale}
  ...
}
\author{ D. Nychka (but really who ever wrote the persp help file) }

\seealso{ persp, drape.plot} 
\examples{
# Mr. R's favorite New  Zealand Volcano!
data( volcano)
M<- nrow( volcano)
N<- ncol( volcano)
x<- seq( 0,1,,M)
y<- seq( 0,1,,N)

persp( x,y,volcano,phi=30, theta=30)-> pm
#point at (.5,.5,.9)
trans3d( .5,.5,180, pm)-> uv
points( uv, pch="o", cex=2, col="red")

# add a contour line 
# 
contourLines( x,y,volcano)[10][[1]]-> out
z<- rep( out$level, length( out$x))
trans3d( out$x, out$y, z, pm)-> uv
lines( uv, col="blue", lwd=2)



}
\keyword{ hplot }% at least one, from doc/KEYWORDS
