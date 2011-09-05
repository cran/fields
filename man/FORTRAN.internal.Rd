\name{fields exported FORTRAN}
 \alias{css}
 \alias{ddfind}
 \alias{dmaket}
 \alias{evlpoly}
 \alias{evlpoly2}
 \alias{igpoly}
 \alias{inpoly}
 \alias{multeb}
 \alias{multrb}
 \alias{radbas}
 \alias{rcss}
\title{
FORTRAN subroutines used in fields functions
}
\description{
These functions implement cubic smoothing splines and also provide some basic 
computations for radial basis functions. All are called using the \code{.FORTRAN}
interface.
}
\details{
For these low level FORTRAN subroutines refer to the R functions for the calling
sequence and to the src subdirectory of the fields pacakage for the source code. 
\describe{
\item{css}{Cubic smoothing spline see \code{sreg} and \code{splint}}
\item{ddfind}{Finds nearest neighbor points within a fixed distance. See \code{fields.rdist.near}}
\item{dmaket}{Creates matrix of all polynomial terms up to fixed order. See \code{fields.mkpoly}}
\item{evlpoly}{evaluates a univariate polynomial. See code{fields.evlpoly}}
\item{evlpoly2}{ evaluates a multivariate polynomial. See code{fields.evlpoly2}}
 \item{inpoly}{Determine which  2-d locations are within a polynomial. see \code{in.poly}}
 \item{igpoly}{Determine which  2-d grid points locations are within a polynomial. see \code{in.poly.grid}}
 \item{multeb}{Multiply an exponential cross covariance matrix by another matrix. See \code{exp.cov}}
  \item{multrb}{Multiply an radial basis function matrix by another matrix. See \code{rad.cov}}
 \item{radbas}{Evaluates radial basis functions. See \code{rdist.R}}
 \item{rcss}{Robust cubic smoothing spline. See \code{qsreg}}
}
}
\keyword{internal}
