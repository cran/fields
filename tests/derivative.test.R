# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

#library( fields, lib.loc="lib.test")

library( fields)
options(echo=FALSE)
test.for.zero.flag<- 1

DD<- cbind( seq(.01,2,,50))
look2<- Wendland(DD, theta=1.0, dimension=2,k=3,derivative=1) 

look1<- (Wendland(DD+1e-5, theta=1.0, dimension=2,k=3)
- Wendland(DD-1e-5, theta=1.0, dimension=2,k=3))/2e-5

test.for.zero( look1, look2,tol=1e-6)

look2<- Wendland(DD, theta=1.5, dimension=2,k=3,derivative=1) 

look1<- (Wendland(DD+1e-5, theta=1.5, dimension=2,k=3)
- Wendland(DD-1e-5, theta=1.5, dimension=2,k=3))/2e-5

test.for.zero( look1, look2,tol=1e-6)

x<- seq( -1,1,,15)

ctest<- rep(0,15)
ctest[8]<- 1

wendland.cov( x,x, k=2, theta=.75)-> look0
look0%*% ctest->look0

wendland.cov( x,x, k=2, theta=.75, C=ctest, derivative=0)-> look

test.for.zero( look0, look)


wendland.cov( x,x, k=2, theta=1.0, C=ctest, derivative=1)-> look

wendland.cov( x+1e-5, x, k=2, theta=1.0, C=ctest)-
wendland.cov( x-1e-5, x, k=2, theta=1.0, C=ctest)-> look2
look2<- look2/2e-5
 
test.for.zero( look, look2,tol=1e-7)


wendland.cov( x,x, k=2, theta=.75, C=ctest, derivative=1)-> look

wendland.cov( x+1e-5, x, k=2, theta=.75, C=ctest)-
wendland.cov( x-1e-5, x, k=2, theta=.75, C=ctest)-> look2
look2<- look2/2e-5

test.for.zero( look, look2,tol=1e-7)



x<- make.surface.grid( list(x=seq( -1,1,,40), y=seq( -1,1,,40)))
#x<- make.surface.grid( list(x=seq( -1,1,,20), y=seq( -1,1,,20)))

y<- (.123*x[,1] + .234*x[,2])
obj<- mKrig( x,y, lambda=0, cov.function="wendland.cov", k=3, theta=.2)
#obj<- mKrig( x,y, lambda=0, m=3, theta=.2)


xp<- make.surface.grid( list(x=seq(-.5,.5,,24),y= seq( -.5,.5,,24)) )
predict( obj, xp, derivative=1)-> outd
test.for.zero( outd[,1],.123)
test.for.zero( outd[,2],.234)

y<- (x[,1]**2 - 2* x[,1]*x[,2] +  x[,2]**2)/2
obj<- mKrig( x,y, lambda=0, cov.function="wendland.cov", k=3, theta=.2)

predict( obj, xp, derivative=1)-> outd


true<- cbind( xp[,1] -  xp[,2], xp[,2]- xp[,1])

rmse<-sqrt(mean((true[,1] - outd[,1])**2)/mean(true[,1]**2))
test.for.zero( rmse,0, tol=1e-2,relative=FALSE)

obj<- mKrig( x,y, lambda=0, cov.function="wendland.cov", k=3, V=diag(c( .2,.2) ))

predict( obj, xp, derivative=1)-> outd
rmse<-sqrt(mean((true[,1] - outd[,1])**2)/mean(true[,1]**2))
test.for.zero( rmse,0, tol=1e-2,relative=FALSE)

obj<- mKrig( x,y, lambda=0, cov.function="wendland.cov", k=3,
             V=diag(c( .3,.3) )%*% matrix( c( .5,.5,.5, -.5), 2,2))

predict( obj, xp, derivative=1)-> outd
rmse<-sqrt(mean((true[,1] - outd[,1])**2)/mean(true[,1]**2))
test.for.zero( rmse,0, tol=1e-2,relative=FALSE)



rmse<-sqrt(mean((true[,2] - outd[,2])**2)/mean(true[,2]**2))
test.for.zero( rmse,0, tol=1e-2,relative=FALSE)

cat("done with dervative tests", fill=TRUE)
options( echo=TRUE)

