\name{world}
\alias{world}
\title{Plot of the world}
\description{
Plots quickly, medium resolution outlines of large land masses and 
bodies of water.  
}
\usage{
world(ylim=c(-90, 90), xlim=c(-180, 180), add=FALSE, asp = 1,
      xlab = "", ylab = "",  xaxt = "n", yaxt = "n", ...)
}
\arguments{
  \item{ylim}{range of latitudes}
  \item{xlim}{range of longitudes}
  \item{add}{logical; if true will add the world map to current plot.}
  \item{asp}{aspect ratio used if \code{add} is false, see
    \code{\link{plot.default}}.}
  \item{xlab,ylab}{labels for x- and y-axis; empty by default.}
  \item{xaxt,yaxt}{axis type for x- and y-axis; empty by default.}
 \item{\dots}{
    These are graphical arguments that are passed to the
    \code{\link{lines}} (and \code{\link{plot}} if \code{add} is false)
    function that draws the outline.}
}
\details{
This function uses the FIELDS dataset world.dat for the coordinates. The 
main advantage of this function is that it is fast and easy to modify. 
 }
\examples{
# Draw map in device color number 2 ( usually red?) 

world( col=2)

## Western Europe (*which* big island is missing ?)
## with a (light color!) coordinate grid:

world(xlim=c(-10,18),ylim=c(36,60), xaxt = "s", yaxt = "s")
grid()

}
\keyword{hplot}
% docclass is function
