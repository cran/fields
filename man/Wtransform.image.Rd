\name{Wtransform.image}
\alias{Wtransform.image}
\title{
  Quadratic W wavelet transform for an image 
}
\description{
Finds the forward or inverse discrete wavelet transform using the
quadratic W basis. 
}
\usage{
Wtransform.image(x, inv=F, transpose=F, cut.min=8)
}
\arguments{
\item{x}{
Matrix to be transformed 
}
\item{inv}{
If true will compute the inverse transform default is false  
}
\item{transpose}{
If true will compute the transpose of transform default is false 
}
\item{cut.min}{
Minimum level of transformation.  cut.min=8 means that the coarsest level
will consist of 64=8X8 scale functions centered on an 8X8 grid.  
}
}
\value{
A matrix the same size as x.   
}
\details{
The wavelet transform can be thought as matrix multiplication 
A \%*\% vec(x)  where vec(x) is the matrix x  stacked by columns.  
The inverse transform is  inv(A) \%*\% vec(x) and transpose is  
 t(A) \%*\% vec( x) 
}
\section{References}{
Nychka,D.  Wikle,C. ,  
Royle, J.A. (1999) Large spatial prediction problems and nonstationary 
random 
fields  
}
\seealso{
plot.coef, WQS, WQSi, Wtransform.D  
}
\examples{
# Wtransform of John Lennon image
data(lennon)
look<- Wtransform.image( lennon)
image.plot( look)
#threshhold 
thr<-  quantile( abs( look), .95)
temp<- look
temp[abs(look)< thr] <- 0
look2<- Wtransform.image( temp, inv=T)
image( look2) # 95 \% compressed image

 
# a diagonal detail basis function 
temp<- matrix(0, nrow=32, ncol=32) 
temp[12,12]<- 1 
look<- Wtransform.image( temp , inv=T)
persp( look) 
}
\keyword{spatial}
% docclass is function
% Converted by Sd2Rd version 1.21.
