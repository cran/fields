\name{The Krig engines and support functions}
\alias{Krig.engine.default}
\alias{Krig.engine.knots}
\alias{Krig.engine.fixed}
\alias{Krig.coef}
\alias{Krig.check.xY}
\alias{Krig.cor.Y}
\alias{Krig.transform.xY}

\title{ Basic linear algebra computations supporting the Krig function. }
\description{
These are internal functions to Krig that compute the basic matrix
decompositions or solve the linear systems needed to evaluate the 
Krig/Tps estimate. Others listed below do some simple housekeeping and 
formatting. Typically they are called from within Krig but can also be
used directly if passed a Krig object list. 
}
\usage{
Krig.engine.default(out, verbose = FALSE)
Krig.engine.knots(out, verbose = FALSE)
Krig.engine.fixed( out, verbose=FALSE, lambda=NA)

Krig.coef(out, lambda = out$lambda, y = NULL, yM = NULL, verbose = FALSE) 
Krig.check.xY(x, Y, weights, na.rm, verbose = FALSE) 
Krig.cor.Y(obj, verbose = FALSE) 
Krig.transform.xY(obj, knots, verbose = FALSE)

}
\arguments{
  \item{out}{ A complete or partial Krig object. If partial it must have all the information accumalted to this calling 
point within the 
Krig function.  }

\item{obj}{ Same as \code{out}. }

\item{verbose}{ If TRUE prints out intermediate results for
debugging.}

\item{lambda}{ Value of smoothing parameter "hard wired" into decompostions. 
Default is NA -- use the value in the out object. }

\item{y} {New y vector for recomputing coefficients}

\item{yM} {New y vector for recomputing coefficients but 
the values have already been collapsed into replicate group means.}

\item{Y} {raw data Y vector}

\item{x} {raw x matrix}

\item{weights} { Raw \code{weights} vector passed to Krig}

\item{na.rm} { NA action logical vlaues passed to Krig}

\item{knots} {Raw \code{knots} matrix  passed to Krig} 

}

\details{

ENGINES:

The engines are the 
code modules that handle the basic linear algebra needed to 
computed the estimated curve or surface coefficients. 
All the engine work on the data that has been reduced to unique 
locations and possibly replicate group means with the weights 
adjusted accordingly. All information needed for the decomposition are 
components in the Krig object passed to these functions. 

 \code{Krig.engine.default} finds the decompositions for a Universal
Kriging estimator. by simultanously diagonalizing the linear system
system for the coefficients of the estimator. The main advantage of this
form is that it is fairly stable numerically, even with ill-conditioned
covariance matrices with lambda > 0. (i.e.  provided there is a "nugget"
or measure measurement error. Also the eigendecomposition allows for
rapid evaluation of the likelihood, GCV and coefficients for new data
vectors under different values of the smoothing parameter, lambda. 


 \code{Krig.engine.knots} finds the decompositions in the case that the
covariance is evaluated at arbitrary locations possibly different than
the data locations (called knots). The intent of these decompostions is
to facilitate the evaluation at different values for lambda.  There will
be computational savings when the number of knots is less than the
number of unique locations. (But the knots are as densely distributed as
the structure in the underlying spatial process.) This function call
fields.diagonalize, a function that computes the matrix and eigenvalues
that simultaneous diagonalize a nonnegative definite and a positive
definite matrix. These decompostions also facilitate multiple
evaluations of the likelihood and GCV functions in estimating a
smoothing parameter and also multiple solutions for different y vectors. 

 \code{Krig.engine.fixed} are specific decomposition based on the Cholesky 
factorization assuming that the smoothing parameter is fixed. This 
is the only case that works in the sparse matrix.
Both knots and the full set of locations can be handled by this case. 
The difference between the "knots" engine above is that only a single value
of lambda is considered in the fixed engine. 


OTHER FUNCTIONS:

\code{Krig.coef} Computes the "c" and "d" coefficients to represent the 
esimated curve. These coefficients are used by the predict functions for 
evaluations. Krig.coef can be used outside of the call to Krig to 
recompute the fit with different Y values and possibly with different
lambda values. If new y values are not passed to this function then the yM
vector in the Krig object is used. The internal function 
\code{Krig.ynew} sorts out the logic of what to do and use based on the 
passed arguments. 

\code{Krig.check.xY} Checks for removes missing values (NAs).

\code{Krig.cor.Y} Standardizes the data vector Y based on a correlation model. 

 \code{Krig.transform.xY} Finds all replcates and collapse to unique
locations and mean response and pooled variances and weights. These are
the xM, yM and weightsM used in the engines. Also scales the x locations
and the knots according to the transformation. 

}

\value{

ENGINES:

The returned value is a list with the matix decompositions and 
other information. These are incorporated into the complete Krig object. 

Common to all engines:

   \item{decomp}{Type of decomposition}
   \item{nt}{ dimension of T matrix}
   \item{np}{ number of knots}

Krig.engine.default:

\item{u}{Transformed data using eigenvectors.}
\item{D}{Eigenvalues}
\item{G}{Reduced and weighted matrix of the eigenvectors}
\item{qr.T}{QR decomposition of fixed regression matrix}
\item{V}{The eigenvectors}


Krig.engine.knots:

  \item{u}{A transformed vector that is based on the data vector.}
  \item{D }{Eigenvalues of decomposition}
  \item{G} { Matrix from diagonalization}
  \item{qr.T} {QR decomposition of the matrix for the fixed component. 
 i.e. sqrt( Wm)\%*\%T}
   \item{pure.ss}{ pure error sums of squares including both the
    variance from replicates and also the sums of squared residuals
    from fitting the full knot model with lambda=0 to the replicate means. }


Krig.engine.fixed:

\item{d} { estimated coeficients for the fixed part of model}
\item{c} { estimated coeficients for the basis functions dervied from the 
                 covariance function.}

Using all data locations

\item{qr.VT}{QR decomposition of the inverse Cholesky factor times the 
T matrix. }
\item{MC}{Cholesky factor}

Using knot locations

\item{qr.Treg}{QR decomposition of regression matrix modified by the 
estimate of the nonparametric ( or spatial) component.} 
\item{lambda.fixed}{ Value of lambda used in the decompositions}

OTHER FUNCTIONS:

\code{Krig.coef}
\item{yM}{Y values as replicate group means}
\item{shat.rep}{Sample standard deviation of replicates}
\item{shat.pure.error} {Same as shat.rep}
\item{pure.ss}{ Pure error sums of squares based on replicates}
\item{c}{ The "c" basis coefficients associated with the covariance
or radial basis functionss.}
\item{d}{ The "d" regression type coefficients that are from the fixed part of the model
or the linear null space.}
\item{u}{When the default decomposition is used the data vecotr transformed by the orthogonal matrices. This facilitates evaluting the GCV function 
at different values of the smoothing parameter.}
              

}
\author{Doug Nychka }
\seealso{  \code{\link{Krig}}, \code{\link{Tps}} }
\examples{

Krig( ozone$x, ozone$y)-> out

Krig.engine.default( out)-> stuff

# compare "stuff" to components in out$matrices

Krig.coef( out)$c
# compare to out$c

Krig.coef( out, yM = ozone$y)$c
# better be the same even though we pass as new data!

}
\keyword{ spatial }
