# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

library( fields)
options( echo=FALSE)

##########################################
  test.for.zero.flag<- 1
  data( ozone2)
  x<-ozone2$lon.lat
  y<- ozone2$y[16,]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]

  N<- length( y)
  beta<- -.2
  lambda <-  1.5
  obj<- LatticeKrig( x,y,NC=16, lambda=lambda, beta=beta, NtrA=20,iseed=122)

  grid.info<- obj$MRinfo$grid.info
  K<- MR.cov( x,x,grid.info, beta=beta)
  tempM<-  K
  diag(tempM) <- (lambda) + diag(tempM)
  Mi<- solve( tempM)
  T.matrix<- cbind( rep(1,N),x) 
  d.coef0 <-  solve( t(T.matrix)%*%Mi%*%T.matrix, t(T.matrix)%*%Mi%*%y)

  test.for.zero( obj$d.coef, d.coef0, tag="d from LatticeKrig and by hand")
#### this is c for standard Kriging equations as done in mKrig
  temp2<- chol( tempM)
  c.coef0 <- forwardsolve(temp2, transpose = TRUE,
                        (y- T.matrix%*%d.coef0), upper.tri = TRUE)
  c.coef0 <- backsolve(temp2, c.coef0)

### find these using mKrig (still standard Kriging) and the lattice covariance function:
  obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="MR.cov",
                                 cov.args=list(grid.info=grid.info, beta=beta),
                                 NtrA=20, iseed=122)
  test.for.zero( obj0$c, c.coef0, tag="c from mKrig and by hand" )
# we also know that for standard Kriging
# residuals = lambda* c.coef0
# use this to check the initial LatticeKrig result
 test.for.zero( obj0$fitted.values, obj$fitted.values)
 test.for.zero( lambda*obj0$c, (y-obj$fitted.values),
               tag="c from mKrig and from residuals of LatticeKrig (this is big!)" )
#
# test more complex covariance model:
#
  alpha<- c(1,.5,.5)
  nlevel<-3
  beta<- c( -.2,-.2,.05)
  obj<- LatticeKrig( x,y,NC=5, lambda=lambda,
                        nlevel=nlevel, alpha=alpha,beta=beta, NtrA=20,iseed=122)
  grid.info<- obj$MRinfo$grid.info
  obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="MR.cov",
                                 cov.args=list(grid.info=grid.info,nlevel=nlevel,
                                 alpha=alpha, beta=beta),
                                 NtrA=20, iseed=122)
  test.for.zero( obj0$fitted.values, obj$fitted.values)
  test.for.zero( obj$d.coef, obj0$d, tag= "d from Lattice Krig and mKrig")
#### compare estimates of trace
  test.for.zero( obj$trA.info, obj0$trA.info, tag="Monte Carlo traces")

# evaluate predicted values
  glist<- fields.x.to.grid( x,10, 10)
  xg<-  make.surface.grid(glist)
  grid.info<- obj$MRinfo$grid.info
  MRinfo<- obj$MRinfo
# first "by hand"
  Tmatrix<- cbind( rep(1,nrow(xg)), xg)
  yhat0<- Tmatrix%*%obj0$d +
            MR.cov( xg,x,grid.info,nlevel=nlevel,alpha=alpha,beta=beta)%*%obj0$c
  
  PHIg<- MRbasis( xg, MRinfo)
  yhat1<- Tmatrix%*%obj$d.coef + PHIg%*%obj$c.coef
  test.for.zero( yhat1, yhat0, tag="predicted values   by hand and  from standard and DRK")
  test.for.zero( yhat1, predict(obj,xg), tag="predicted values from LatticeKrig and by hand"  )
  test.for.zero( predict(obj,xg), predict(obj0,xg),
                           tag="predicted values LatticeKrig and mKrig")

# tests for computing the determinant and quad form
  test.for.zero.flag<- 1
  alpha<- c(1,.5,.5)
  nlevel<-3
  beta<- c( -.2,-.2,.05)

  lnDet<- function(A){
  sum( log( eigen( A, symmetric=TRUE)$values))}

  data( ozone2)
  x<-ozone2$lon.lat[1:20,]
  y<- ozone2$y[16,1:20]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
#x<- transformx(x, "range")
  N<- length( y)
  lambda <- .8
# a micro sized lattice so determinant is not too big or small
  obj<- LatticeKrig( x,y,NC=5, lambda=lambda,nlevel=nlevel,alpha=alpha,
                              beta=beta, NtrA=5,iseed=122)
  MRinfo<- obj$MRinfo
  grid.info<- MRinfo$grid.info
  PHI<- MRbasis( x,MRinfo, spam.format=TRUE)
  Q <- MRprecision(MRinfo, spam.format = TRUE, beta = beta, alpha=alpha)
# coerce to full matrix
  Q<- spam2full(Q)
  Mtest<- PHI%*% (solve( Q)) %*% t( PHI) + diag(lambda, N)
  temp<- t(PHI)%*%PHI + lambda*Q
  A<- Q*lambda
  B1<-  PHI%*% (solve( A)) %*% t( PHI) + diag(1, N)
  B2<-  t(PHI)%*%PHI + A
# the bullet proof application of identity
  test.for.zero(lnDet( B1),lnDet( B2)- lnDet(A))
# now adjusting for lambda factor 
  Mtest<- PHI%*% (solve( Q)) %*% t( PHI) + diag(lambda, N)
  test.for.zero( lambda*B1, Mtest)
  test.for.zero(lnDet( Mtest), lnDet(B2) - lnDet(lambda*Q) + N*log(lambda) )
  test.for.zero(lnDet( Mtest), lnDet(B2) - lnDet(Q) + (-MRinfo$Ntotal + N)*log(lambda) )

# find log determinant of temp using cholesky factors
# applying det identity
   temp<- t(PHI)%*%PHI + lambda*Q
   chol( temp)-> Mc

   lnDetReg <- 2 * sum(log(diag(Mc)))
   lnDetQ<-  2* sum( log( diag( chol(Q))))
   lnDetCov<- lnDetReg - lnDetQ + (-MRinfo$Ntotal + N)*log(lambda)
   test.for.zero( lnDetCov, lnDet( Mtest))
   test.for.zero( obj$lnDetCov, lnDet( Mtest), tag="LatticeKrig and direct test of lnDetCov")
#
# now check these formulas as implemented in LatticeKrig
    obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="MR.cov",
                                 cov.args=list(grid.info=grid.info,nlevel=nlevel,
                                 alpha=alpha, beta=beta),
                                 NtrA=20, iseed=122)
 
 test.for.zero( obj$lnDetCov,obj0$lnDetCov, tag= "lnDetCov for mKrig and LatticeKrig")
 test.for.zero( obj$quad.form,  obj0$quad.form, tag= "quadratic forms for rho hat")
 test.for.zero(  obj0$lnProfileLike, obj$lnProfileLike,
                                tag="Profile Likelihood concentrated on lambda" )

# repeat tests for weighted measurement errors.
# recopy data to make reading easier
  data( ozone2)
  x<-ozone2$lon.lat[1:20,]
  y<- ozone2$y[16,1:20]
  good <-  !is.na( y)
  x<- x[good,]
  y<- y[good]
#x<- transformx(x, "range")
  N<- length( y)
 alpha<- c(1,.5,.5)
  nlevel<-3
  beta<- c( -.2,-.2,.05)
  lambda <- .8
  N<- length(y)
  set.seed(243)
  weights<- runif(N)*10
  test.for.zero.flag<- 1
  obj<- LatticeKrig( x,y,weights,NC=15,
                    lambda=lambda,alpha=alpha,nlevel=nlevel, beta=beta, NtrA=5,iseed=122)
        grid.info<- obj$MRinfo$grid.info
# compare mKrig and Krig with weights and LatticeKrig
  obj0<- mKrig( x,y,weights, lambda=lambda, m=2, cov.function="MR.cov",
                                 cov.args=list(grid.info=grid.info,nlevel=nlevel,
                                 alpha=alpha, beta=beta),
                                 NtrA=20, iseed=122)
 
  obj1<- Krig( x,y,weights=weights, lambda=lambda,GCV=TRUE, m=2,
               cov.function="MR.cov", cov.args=list(grid.info=grid.info,nlevel=nlevel,
               alpha=alpha, beta=beta))
            
 test.for.zero( obj0$fitted.values, obj1$fitted.values)
 test.for.zero( predict(obj0), predict(obj1), tag="predicted  values mKrig/Krig  w/weights")
 test.for.zero( obj0$rhohat, obj1$rhohat,tag="compare rhohat for mKrig and Krig with weights")

############ now tests for LatticeKrig

 test.for.zero( obj$fitted.values, obj0$fitted.values)
############# tests using setup option

 obj2<- LatticeKrig( x,weights=weights,NC=15, lambda=lambda,alpha=alpha,
                    nlevel=nlevel,beta=beta, setup=TRUE)
 look<- LatticeKrig.coef( obj2$Mc, obj2$wPHI, obj2$wT.matrix, y*sqrt(weights), lambda)

 test.for.zero( look$c.coef, obj$c.coef, tag="test of LatticeKrig.coef c")
 test.for.zero( look$d.coef, obj$d.coef, tag="test of LatticeKrig.coef d")

 Q<- MRprecision(obj2$MRinfo, spam.format=TRUE,alpha=alpha, beta=beta)
 look2<-LatticeKrig.lnPlike(obj2$Mc,Q,y, 
                             lambda,obj$residuals, weights)
 test.for.zero( look2$lnProfileLike, obj$lnProfileLike)

# all done!
 cat("Done testing LatticeKrig",fill=TRUE)
 options( echo=FALSE)


# SUPPLEMENT: commented out sanity checks for  weighted/unweighted versions of mKrig and Krig
#hold0<-Krig ( x,y,weights=weights,method="user",GCV=TRUE,lambda=1e-3,
#             cov.function="Exp.simple.cov", cov.args=list( theta=300) )
#hold1<-mKrig(x,y,weights, lambda=1e-3,cov.function="Exp.simple.cov", cov.args=list( theta=300))
#test.for.zero( predict(hold0), predict(hold1))

# obj1<- Krig( x,y,  lambda=lambda,GCV=FALSE, m=2, cov.function="LatticeKrig.cov",
#               cov.args=list( center.grid.list=obj$center.grid.list,beta=beta))

# obj0<- mKrig( x,y, lambda=lambda, m=2, cov.function="LatticeKrig.cov",
#               cov.args=list( center.grid.list=obj$center.grid.list,beta=beta),
#               NtrA=5, iseed=122)
# test.for.zero( predict(obj0), predict(obj1))
# test.for.zero( obj0$shat.MLE, obj0$shat.MLE)
#








