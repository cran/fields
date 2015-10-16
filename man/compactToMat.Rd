\name{compactToMat}
\alias{compactToMat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Convert Matrix from Compact Vector to Standard Form
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
\code{compactToMat} transforms a matrix from compact, vector form to
  a standard matrix.  Only symmetric matrices can be stored in this
  form, since a compact matrix is stored as a vector  with elements 
  representing the upper triangle of the matrix.  This function assumes
  the vector does not contain diagonal elements of the matrix.
  
  An example of a matrix stored in compact form is any matrix 
  generated from the \code{rdist} function with \code{compact=TRUE}.
}
\usage{
compactToMat(compactMat, diagVal=0, lower.tri=FALSE, upper.tri=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{compactMat}{
%%     ~~Describe \code{compactMat} here~~
A symmetric matrix stored as a vector containing elements for the lower-triangular
  portion of the true matrix (and none of the diagonal elements), as returned by 
  \code{rdist} with \code{compact=TRUE}.
}
  \item{diagVal}{
%%     ~~Describe \code{diagVal} here~~
A number to put in the diagonal entries of the output matrix.
}
  \item{lower.tri}{
%%     ~~Describe \code{diagVal} here~~
Whether or not to fill in the upper triangle of the output matrix
}
  \item{upper.tri}{
%%     ~~Describe \code{diagVal} here~~
Whether or not to fill in the lower triangle of the output matrix
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
The standard form matrix represented by the input compact matrix
}
\author{
%%  ~~who you are~~
John Paige
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{rdist}}, \code{link{dist}}
}
\examples{
################
#Calculate distance matrix from compact form:
################

#make a distance matrix
distOut = rdist(1:5, compact=TRUE)
print(distOut)

#note that distOut is in compact form:
print(c(distOut))

#convert to standard matrix form:
distMat = compactToMat(distOut)

################
#fast computation of covariance matrix:
################

#generate 5 random points on [0,1]x[0,1] square
x = matrix(runif(10), nrow=5)

#get compact distance matrix
distOut = rdist(x, compact=TRUE)

#evaluate Exponential covariance with range=1.  Note that
#Covariance function is only evaluated over upper triangle
#so time is saved.
diagVal = Exponential(0, range=1)
compactCovMat = Exponential(distOut, range=1)
upperCovMat = compactToMat(compactCovMat, diagVal)
lowerCovMat = compactToMat(compactCovMat, diagVal, lower.tri=TRUE, upper.tri=FALSE)
fullCovMat = compactToMat(compactCovMat, diagVal, lower.tri=TRUE, upper.tri=TRUE)
compactCovMat
lowerCovMat
upperCovMat
fullCovMat
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ compact }
\keyword{ matrix }% __ONLY ONE__ keyword per line
