##################################
######### testing functions
##################################

MaternGLS.test<- function(x,y, smoothness=1.5,init= log( c(1.0,.2,.1))){
# some simulations within fields to
# study the variability in estimates of the covariance parameters.
    N<- length( y)
    Tmatrix <- fields.mkpoly(x,2)
    qr.T <- qr(Tmatrix)
    Q2 <- qr.yq2(qr.T, diag( 1,N))
    nu<- smoothness
    
loglmvn <- function(pars,nu,y,x,d){
	lrho=pars[1]; ltheta=pars[2]; lsig2=pars[3]
     #  print( pars)
	N <-length(y)
       	M <- ( exp(lrho)*Matern(d,range=exp(ltheta),smoothness=nu) + 
                          exp(lsig2)*diag(N) )
        
X<- fields.mkpoly( x, 2)
Mi<- solve( M)
betahat<-  solve(t(X)%*%Mi%*%X)%*% t(X)%*% Mi%*% y
res<- y - X%*%betahat
chol(M)-> cM
lLike<-  (N/2)*log(2*pi) - (1/2)* (2*sum( log( diag(cM)))) - (1/2)* t(res)%*% Mi %*% res

        ycept<- -lLike
        
	if((abs(lrho)>20)| (abs(ltheta)>10)|  (abs(lsig2)>40)) {
           return( ycept + 1000*sum(abs(abs(pars)-100))) }
        else{
            return( ycept)}
  }
    
    d<- rdist( x,x)
    temp<- optim(init,loglmvn,method="L-BFGS-B",nu=nu,y=y,x=x,d=d)
    out<- exp(temp$par)

return( list(  smoothness=smoothness, pars= out, optim=temp))
}

 

MaternGLSProfile.test<- function(x,y, smoothness=1.5,init= log( c( .05,1.0))){
# some simulations within fields to
# study the variability in estimates of the covariance parameters.
    N<- length( y)
   
    nu<- smoothness
    
loglmvn <- function(pars,nu,y,x,d){
	llam=pars[1]; ltheta=pars[2]
      # print( pars)
	N <-length(y)
        theta<- exp( ltheta)
        lambda<- exp(llam)
         lLike<- mKrig( x,y, theta=theta,Covariance="Matern", smoothness=nu, lambda=lambda)$lnProfileLike
        ycept<- -lLike
	if((abs(llam)>20)| (abs(ltheta)>10) ) {
           return( ycept + 1000*sum(abs(abs(pars)-100))) }
        else{
            return( ycept)}
  }
    
    d<- rdist( x,x)
    temp<- optim(init,loglmvn,method="L-BFGS-B",nu=nu,y=y,x=x,d=d)
    out<- exp(temp$par)
    rho.MLE<-mKrig( x,y, theta=out[2],Covariance="Matern", smoothness=nu, lambda=out[1])$rho.MLE
    sigma2.MLE<- out[1]*rho.MLE
return( list(  smoothness=smoothness,pars=c(rho.MLE, out[2], sigma2.MLE),
                  optim=temp))
}

MaternQR.test<- function(x,y, smoothness=1.5,init= log( c(1.0,.2,.1))){
# some simulations within fields to
# study the variability in estimates of the covariance parameters.
   nu<- smoothness
    
loglmvn <- function(pars,nu,x,y){
        N<- length( y)
        Tmatrix <- fields.mkpoly(x,2)
        qr.T <- qr(Tmatrix)
        Q2 <- qr.yq2(qr.T, diag( 1,N))
        ys<- t(Q2)%*%y
        N2<- length( ys)
	lrho=pars[1]; ltheta=pars[2]; lsig2=pars[3]
        d<- rdist( x,x)
        
       	A <-  ( Matern(d,scale=exp(lrho),range=exp(ltheta),smoothness=nu) + exp(lsig2)*diag(N) )
        A <- t(Q2)%*%A%*%Q2
        A<- chol( A)
        w = backsolve(A,ys,transpose=TRUE)
	ycept<- (N2/2)*log(2*pi) + sum(log(diag(A))) + (1/2)*t(w)%*%w
        
	if((abs(lrho)>100)| (abs(ltheta)>100)|  (abs(ltheta)>100)) {
           return( ycept + 1000*sum(abs(abs(pars)-100))) }
        else{
            return( ycept)}
  }
    
   
    temp<- optim(init,loglmvn,method="L-BFGS-B",nu=nu,x=x,y=y)
    out<- exp(temp$par)
    llike<- loglmvn( temp$par, nu,x,y)
return( list(  smoothness=smoothness,pars= out, llike=llike, optim=temp))
}


MaternQRProfile.test<- function(x,y, smoothness=1.5,init= log( c(1.0)) ){
# some simulations within fields to
# study the variability in estimates of the covariance parameters.
  nu<- smoothness
loglmvn <- function(pars,nu,x,y){
	ltheta=pars[1]
       # print( exp(ltheta))
        ycept<- Krig( x,y, Covariance="Matern", theta= exp( ltheta),
                     smoothness=nu, method="REML")$lambda.est[6,5]
#        print( c(exp(ltheta),ycept))
	if((abs(ltheta)>100)){
            return( ycept + 1000*sum(abs(abs(pars)-100))) }
        else{
            return( ycept)}
  }
    
    temp<- optim(init,loglmvn,method="L-BFGS-B",nu=nu,x=x,y=y)
    theta.est<- exp(temp$par[1])
    out2<-Krig( x,y, Covariance="Matern", theta= theta.est,smoothness=nu, method="REML")
# MLE based on reduced degrees of freedom:

    offset<- (out2$N/( out2$N-3))
    
    out3<- c( out2$rho.MLE*offset, theta.est, out2$shat.MLE^2*offset)  

return( list(obj=out2,  smoothness=smoothness,pars= out3, trA= out2$eff.df, optim=temp))
}

# this function has correct formula for REML likelihood
  REML.test<-function (x,y,rho,sigma2, theta, nu=1.5){
  Tmatrix <- fields.mkpoly(x,2)
  qr.T <- qr(Tmatrix)
  N<- length( y)
  Q2 <- qr.yq2(qr.T, diag( 1,N))
  ys<- t(Q2)%*%y
  N2 <- length( ys)
  A <-  (rho* Matern(rdist(x,x),range=theta,smoothness=nu) + sigma2*diag(1,N) )
  A <- t(Q2)%*%A%*%Q2
  Ac<- chol( A)
  w <- backsolve(Ac,ys,transpose=TRUE)
  REML.like<- (N2/2)*log(2*pi) + (1/2)*2*sum(log(diag(Ac))) + (1/2)*t(w)%*%w
  REML.like<- -1* REML.like
  ccoef<-  rho* Q2%*%solve( A)%*%ys
  return( list(REML.like= REML.like, A=A, ccoef=ccoef,
                    quad.form=t(w)%*%w, rhohat= (t(w)%*%w/N2)*rho, det=2*sum(log(diag(Ac))), N2=N2 ) ) 
}



MLE.Matern<- function(x,y,smoothness, theta.grid=NULL,
                       ngrid=20, verbose=FALSE,m=2,niter=25, tol=1e-5,... ){
# remove missing values and print out a warning  
  bad<- is.na( y)
  if( sum(bad)>0){
    cat("removed ",sum(bad)," NAs", fill=TRUE)
    x<- x[!bad,]
    y<- y[!bad] }

  
   objective.fn<- 
     function( ltheta, info){
     minus.lPLike<- Krig( info$x,info$y,Covariance="Matern",
       smoothness= info$smoothness, theta=exp(ltheta), 
       method="REML", nstep.cv=80, give.warnings=FALSE,m=m,...)$lambda.est[6,5]
   return(minus.lPLike)} 

# list to pass to the objective function
 
    info<- list( x=x, y=y, smoothness= smoothness)
 
#
# if grid for ranges is missing use some quantiles of pairwise distances among data.
#  this will only work if the likelihood at endpoints is smaller than middle.
# (i.e. convex)
  
  if( is.null(theta.grid)){
    theta.range<- quantile(  rdist( x,x), c( .03,.97))
    theta.grid <-  seq( theta.range[1], theta.range[2],,ngrid)}
   if( length( theta.grid)==2){
      theta.grid <-  seq( theta.grid[1], theta.grid[2],,ngrid)}
 
  ngrid<- length( theta.grid)
  sighat<- rhohat<- trA<-theta<- rep( NA,ngrid)
  minus.REML<- rep( NA, ngrid)

# grid search
    for( j in 1:ngrid){
   minus.REML[j] <- objective.fn( log(theta.grid[j]), info)       

        }
     IMIN<- ( 1: ngrid)[ min(minus.REML)==minus.REML]

# best point for theta from grid search
  if( IMIN ==1 |IMIN ==ngrid){
    cat("REML at end of search interval:", fill=TRUE)
     temp<- cbind(theta.grid,-minus.REML)
    dimnames(temp)<- list( NULL, c("theta", "logProfileLike"))
    return( list(smoothness=smoothness, pars= rep( NA,3) ,
      REML=NA, trA=NA,REML.grid= temp))}
                                 
   lstart<- log( theta.grid)[ IMIN + c( -1,0,1)]
# golden  section search -- this assumes convex minus log likelihood
# note that search is in log scale.    
  out<- golden.section.search( lstart[1], lstart[2], lstart[3], f=objective.fn, f.extra=info,
                              niter=niter, tol=tol)$x
  theta.MLE<- exp(out)

# one final call to Krig with the theta.MLE value to recover MLEs for rho and sigma
   
  hold<- Krig( info$x,info$y,Covariance="Matern",
             smoothness= info$smoothness, theta=theta.MLE,
             method="REML")
 
  sig.MLE<- hold$shat.MLE
  rhohat<- hold$rhohat
  trA<- hold$lambda.est[6,2]
  REML<- hold$lambda.est[6,5]  # actually minus log REML 
 
  out<- c( rhohat, theta.MLE,sig.MLE)
  names(out)<- c( "rho","theta", "sigma")
  temp<- cbind(theta.grid,-minus.REML)
  dimnames(temp)<- list( NULL, c("theta", "logProfileLike"))
return( list(smoothness=smoothness, pars= out,REML=-REML, trA=trA,REML.grid= temp) )
       

}
