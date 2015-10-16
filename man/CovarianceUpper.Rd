\name{CovarianceUpper}
\alias{ExponentialUpper}
\title{
Evaluate covariance over upper triangle of distance matrix
}
\description{
Evaluates the covariance over the upper triangle of a distance matrix 
  rather than over the entire matrix to reduce computation time.
}
\usage{
ExponentialUpper(distMat, range = 1, alpha = 1/range, phi = 1)
}
\arguments{
  \item{distMat}{
The distance matrix to evaluate the covariance over.
}
  \item{range}{
Range parameter default is one. Note that the scale can also
  be specified through the "theta" scaling argument used in
  fields covariance functions)
}
  \item{alpha}{
1/range
}
  \item{phi}{
Marginal variance.
}
}
\value{
The covariance matrix, where only the upper triangle is calculated.
}
\author{
John Paige
}
\seealso{
\code{\link[fields]{Exponential}}
}
\keyword{ covariance }
