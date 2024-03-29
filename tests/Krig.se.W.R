#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2022 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2
##END HEADER
##END HEADER


suppressMessages(library(fields))
# tests of predictSE using 
# off diag weight matrix for obs (W)

options( echo=FALSE)

test.for.zero.flag<- 1

# a nasty example with off diagonal weights.

set.seed(123)

N<- 50 
x<- matrix( runif( N*2), N,2)
y<- rnorm( N)*.2 + 2*x[,1]**2 +  x[,2]**2 


weights<- runif(N)*10
x0<- cbind( c(.1,.2,.6,.65,.8), c(.05,.5,.73,.9,.95))



temp.wght<- function(x, alpha=.3){
  Exp.cov( x, aRange=.1) }

Krig( x,y, cov.function=Exp.cov,weights=weights,
     wght.function= "temp.wght")-> out
Krig( x,y, cov.function=Exp.cov,weights=weights,W= out$W)-> out2


# direct calculation test for A matrix
#

Krig.Amatrix( out, x=x0)-> A
test.for.zero( A%*%y, predict( out, x0),tag="Amatrix vs. predict")

# now find se.

W2<-out$W2
W<- out$W

Sigma<- out$sigmahat*Exp.cov( out$x,out$x)
temp0<- out$sigmahat*(Exp.cov( x0, x0))
S1<- out$sigmahat*Exp.cov( out$x, x0)

#yhat= Ay
#var( f0 - yhat)=    var( f0) - 2 cov( f0,yhat)+  cov( yhat)

Sigma.obs<-  Krig.make.Wi( out)$Wi
Sigma.obs <- Sigma.obs* (out$tauHat.MLE**2) 

temp1<-  A%*%S1
temp2<- A%*% ( Sigma.obs+ Sigma)%*% t(A)
look<- temp0 - t(temp1) - temp1 +  temp2


#compare to 
# diagonal elements

test<- predictSE( out, x= x0) 
test.for.zero( sqrt(diag(  look)), test,tag="Marginal predictSE")


test<- predictSE( out, x= x0, cov=TRUE)
test2<- predictSE( out2, x= x0, cov=TRUE)
test.for.zero( look, test,tag="Full covariance predictSE")
test.for.zero( look, test2,tag="Full covariance predictSE explicit W")

cat( "all done", fill=TRUE)
options( echo=TRUE)
