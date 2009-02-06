# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


# test of sreg and related functions

library( fields)
options(echo=FALSE)
test.for.zero.flag<- 1

data(ozone2)
x<- ozone2$lon.lat
y<- ozone2$y[16,]

out<- Tps( x,y)

out2<- Krig( x,y, Covariance="RadialBasis", 
           M=2, dimension=2, scale.type="range")








