\name{RMprecip}
\alias{RMprecip}
\title{
  Standardized monthly precipitation for August 1963 inteh Rocky Mountain 
Region
}
\description{
This is a useful spatial data set of moderate size consisting of 800 
locations The list has 
components:
}
\usage{
data(RMprecip)
}
\format{
This data set is a list containing the following components: 
}
\describe{
\item{x}{
Longitude-latitude postion of monitoring stations. 
}
\item{y}{
Standardized precipitaton. The square root of precip has 
been standardized by subtracting the 30 year ( 1960-1989) mean and 
dividing 
by the standard deviation. This is a simple way to transform to a process  
that closer to a stationary Gaussian field.
}
\item{x.cv}{
Additional locations for cross-validation. The 
locations weere generated as a space-filling type design. 
}
\item{y.cv}{
values for cross validation
}
}
\keyword{datasets}
% docclass is data
% Converted by Sd2Rd version 1.21.
