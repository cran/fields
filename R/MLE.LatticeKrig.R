# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
MLE.LatticeKrig<- function(x, y, weights=rep( 1, length(y)),NC, grid.info=NULL,
                           alpha=1.0, beta,lambda, eflag,
                          iseed=123,NtrA=20,
                           verbose=FALSE){
  N<- length(y)
# dummy call to LatticeKrig to setup
  obj<- LatticeKrig(  x,y,NC=NC,grid.info=grid.info, setup=TRUE)
  MRinfo<- obj$MRinfo
# Spatial drift matrix -- assumed to be linear in coordinates. (m=2)
  pg<- as.matrix(expand.grid(alpha,beta, lambda))
  XtX<-  t(obj$wPHI)%*%obj$wPHI
  Q<-  MRprecision(MRinfo,alpha=alpha,beta=beta)
  wy<- weights*y
  NG<- nrow(pg)
# initialize arrays for loop  
  lnProfileLike<-rho.MLE<-shat.MLE<-trA<-SEtrA<-GCV<-rep( NA,NG)
  for (  k in 1: NG){
    ltemp<- pg[k,3]
    btemp<- pg[k,2]
    atemp<- pg[k,1]
    if( verbose){
      cat( "alpha, beta lambda", atemp,  btemp, ltemp, ":", fill=TRUE)}
    Q<-  MRprecision(MRinfo,alpha=atemp,beta=btemp)
    temp<- XtX + ltemp*(Q)
    Mc<-update.spam.chol.NgPeyton(obj$Mc, temp)
   
    out1<- LatticeKrig.coef( Mc, obj$wPHI, obj$wT.matrix, wy, ltemp)
    fitted.values<- (obj$wT.matrix%*%out1$d.coef + obj$wPHI%*%out1$c.coef)/sqrt(weights)
    residuals<- y- fitted.values

  out2<-LatticeKrig.lnPlike(Mc, Q, y,ltemp,residuals, weights)
  
  rho.MLE[k] <- out2$rho.MLE
  shat.MLE[k] <-out2$shat.MLE
  lnProfileLike[k] <-out2$lnProfileLike
    
  if(verbose){
    cat( "rho shat, lnProfileLike", rho.MLE[k],  shat.MLE[k], lnProfileLike[k], fill=TRUE)}
 if(!is.null(NtrA)){
   set.seed(iseed)
   wEy<- matrix( rnorm( NtrA*N), N,NtrA)*sqrt(weights)
   out3<- LatticeKrig.coef( Mc, obj$wPHI, obj$wT.matrix, wEy, ltemp)
   wEyhat<-(obj$wT.matrix%*%out3$d.coef + obj$wPHI%*%out3$c.coef) 
   trA.info <- t((wEy *wEyhat)/weights)%*%rep(1,N)
   trA[k] <- mean(trA.info)
   SEtrA[k]<- sd(trA.info)/sqrt(NtrA)
   GCV[k]<- (sum(residuals**2)/N)/ (( 1- trA[k]/N)**2)
   if(verbose){
     cat("trA", trA[k], fill=TRUE)}
  }
    

} 
  out<- cbind(pg,trA,SEtrA, shat.MLE, rho.MLE, lnProfileLike, GCV)
  dimnames(out)<- list( NULL, c("alpha", "beta", "lambda", "trA", "SEtrA", "sigmaMLE",
                                "rhoMLE","lnProfileLike", "GCV"))
  return( out)
}  
