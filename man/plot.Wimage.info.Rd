\name{plot.Wimage.info}
\alias{plot.Wimage.info}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{plot to check 2-d multiresolution indexing.}
\description{
Plots in image format the differnt blocks of coefficents asscoatied with 
a 2-d multiresolution. This function can be used to check the Wtransfomr indexing 
functions and also as an introduction to how the coefficients are organized. 
}
\usage{
plot.Wimage.info(m, n, cut.min)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{Number of rows in image}
  \item{n}{Number of columns of image }
  \item{cut.min}{Smallest number of smooth basis functions. Coarsest resolution will 
have atleast cut.min X cut.min basis functions with support on a regular grid. }
}
\details{
This function was used to check the (compicated) indexing functions work. 
But it also might a useful graphics to desribe how the different levels of wavelet 
coefficients are pack into the image format. 
}
\author{Doug Nychka}
\seealso{ 
\code{\link{Wtransform.image}}, 
\code{\link{W.info}}, 
\code{\link{Wimage.info}}, 
}
\examples{
#
# 64X 128 image coarsest level has 8X16 smooth basis funcitons. 
#
# NOTE as a matrix the image plot is upside down! Rows are along X-axis and 
# columns on Y
# e.g. the (1,1) element is the lower left corner. 
#

plot.Wimage.info( 64, 128, cut.min=8)

}
\keyword{spatial}% at least one, from doc/KEYWORDS
