% fields, Tools for spatial data
% Copyright 2004-2007, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{Wendland}
\alias{Wendland}
\alias{wendland2.2}
\alias{wendland.coef}
\title{Wendland family of covariance functions}

\description{
 Computes the compactly supported, stationatry Wendland covariance
function as a function ofdistance. This family is useful for creating
sparse covariance matrices. 

}

\usage{
Wendland(d, theta = 1, dimension, k,derivative=0, phi=1)

wendland2.2(d, theta=1)

wendland.coef(d,k)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{d}{Distances between locations. Or for wendland.coef the dimension of
the locations.}
\item{theta}{Scale for distances. This is the same as the
range parameter.}
\item{dimension}{Dimension of the locations}
\item{k}{Order of covariance function.}
\item{derivative}{Indicates derivative of covariance function}
\item{phi}{Scale factor to multiply the function. Equivalent to the
marginal variance or sill if viewed as a covariance function.}
}
\details{
 This is the basic function applied to distances and called by the
\code{wendland.cov} function. It can also be used as the Covariance or 
Taper specifications in the more general 
stationary.cov and station.taper.cov functions.
The Wendland covariance function is a
polynomial on [0,theta] and zero beyond theta. 
The parameter \code{k} detemines the smoothness of the covariance at zero. 
The polynomial
coefficients are computed recursively based on the values of \code{k}
and \code{dimension} in the function \code{wendland.coef}. The
polynomial is evaluated using \code{fields.evlpoly}. 

A specific example of the Wendland family is \code{wendland2.2} and this 
is included mainly for testing. 
}

\value{
A vector of the covariances or its derivative. 

}
\author{Doug Nychka}
\seealso{ wendland.cov, stationary.taper.cov}
\examples{

DD<- seq( 0,1.5,, 200)

y<- Wendland( DD, k=2, dimension=2)

plot( DD, y, type="l")

# should agree with 

y.test<- wendland2.2( DD)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
