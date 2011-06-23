# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


LatticeKrig<- function( x,y=NULL, weights = rep(1, nrow(x)),Z=NULL,NC,
                        lambda=1.0,
                        grid.info=NULL,
                        alpha=1.0, beta=-.2, nlevel=1, eflag=FALSE,                                     iseed=123,NtrA=20, setup=FALSE, verbose=FALSE){

# make sure locations are a matrix and get the rows  
  x<- as.matrix(x)
  N<- nrow(x)
  if (any(duplicated(cat.matrix(x)))) 
        stop("locations are not unique see help(LatticeKrig) ")
# makes sure there are no missing values
  if(!is.null(y)){
    if (any(is.na(y))) 
      stop("Missing values in y should be removed")}
# make sure covariate is a matrix  
   if(!is.null(Z)){
      Z<- as.matrix(Z)}
# if center range is missing use the locations
  if( is.null(grid.info)){
    grid.info<- list( xmin= min( x[,1]), xmax= max( x[,1]),
                      ymin= min( x[,2]), ymax= max( x[,2]))
# spacing for grid
    d1<- grid.info$xmax- grid.info$xmin
    d2<- grid.info$ymax- grid.info$ymin
# actual number of grid points is determined by the spacing delta
# determine delta so that centers are equally
# spaced in both axes and NC is the maximum number of grid points
# along the larger range.     
    grid.info$delta<- max(c(d1,d2))/(NC-1)
  }
# hard wire overlap to 2.5  
  overlap<- 2.5
# hardwire INFLATE to 1.5  
  INFLATE<-1.5
  MRinfo<- MRbasis.setup( grid.info, nlevel=nlevel, eflag=eflag,
                            overlap=overlap,INFLATE=INFLATE)
  if(verbose){
    print(MRinfo)}

#  center.grid.list<-  list( x= seq( x1[1], x1[2],delta ), y= seq( x2[1], x2[2],delta))
#  N1<- length( center.grid.list$x)
#  N2<- length( center.grid.list$y)
#  np<- N1*N2
   np<-  MRinfo$Ntotal
# create centers there should be N1*N2 of these
#  centers<- make.surface.grid( center.grid.list)
# delta used to determine support of wendland basis functions
   N1<- MRinfo$N1
   N2<- MRinfo$N2
#  basis.delta<- overlap*delta
# weighted observation vector  
  wy<- sqrt(weights)*y
# Spatial drift matrix -- assumed to be linear in coordinates. (m=2)
# and includes Z covariate if not NULL  
  wT.matrix<- sqrt(weights)*cbind( rep(1, N), x,Z)
  nt<- ncol(wT.matrix)
  nZ<- ifelse (is.null(Z), 0, ncol(Z))   
  ind.drift <- c(rep(TRUE, (nt-nZ)), rep( FALSE, nZ))
# Matrix of sum( N1*N2) basis function (columns) evaluated at the N locations (rows)
# and multiplied by square root of diagonal weight matrix  
# this can be a large matrix if not encoded in sparse format.
  wS<-  diag.spam(sqrt( weights))
  wPHI<-  wS%*% MRbasis(x,MRinfo) 
  if(verbose){
    print( dim( wPHI))}
# square root of precision matrix of the lattice process
#   solve(t(H)%*%H) is proportional to the covariance matrix of the Markov Random Field
  Q<-  MRprecision(MRinfo, spam.format=TRUE,alpha=alpha,beta=beta)
  
# variational matrix used to find coefficients of fitted surface and evaluate
  
############################################################################################
# this is the regularized regression matrix that is the key to the entire algorithm:
########################################################################################
  temp<- t(wPHI)%*%wPHI + lambda*(Q)
  nzero <- length(temp@entries)
  if(verbose){
    print(nzero)}
#  
# S-M-W identity can be used to evaluate the data covariance matrix:
#   M =  PHI%*% solve(t(H)%*%H) %*% t( PHI) + diag( lambda, N)
# i.e. because temp is sparse, sparse computations can be used to indirectly solve
# linear systems based on M
############################################################################################  
# find Cholesky square root of this matrix
############################################################################################  
#  This is where the heavy lifting happens! Note that temp is sparse format so
#  implied by overloading is a sparse cholesky decomposition. 
#  if this function has been coded efficiently this step should dominate
#  all other computations. 
  chol( temp)-> Mc
# If this is just to set up calculations return the intermeidate results
# this list is used in MLE.LatticeKrig
  if( setup){
     return( list(MRinfo=MRinfo, overlap=overlap, Mc=Mc, np=np,nonzero.entries=nzero,
                  wPHI=wPHI, wT.matrix=wT.matrix, weights=weights, nZ=nZ,
                  ind.drift=ind.drift))
   }
  
  out1<- LatticeKrig.coef( Mc, wPHI, wT.matrix, wy, lambda)
  if( verbose){
    print( out1$d.coef)}
  fitted.values<- (wT.matrix%*%out1$d.coef + wPHI%*%out1$c.coef)/sqrt(weights)
  residuals<- y- fitted.values

  out2<-LatticeKrig.lnPlike(Mc, Q, y,lambda,residuals, weights)
# save seed if random number generation happening outside LatticeKrig
  if(exists(".Random.seed",1)){
    save.seed<- .Random.seed}
  else{
    save.seed <- NA} 
    if( !is.na(iseed)){
      set.seed(iseed)}
# generate N(0,1)  
  wEy<- matrix( rnorm( NtrA*N), N,NtrA)*sqrt(weights)
# restore seed   
   if( !is.na(iseed) & !is.na(save.seed[1])){
         assign(".Random.seed",save.seed, pos=1)}
#  
  out3<- LatticeKrig.coef( Mc, wPHI, wT.matrix, wEy, lambda)
  wEyhat<-(wT.matrix%*%out3$d.coef + wPHI%*%out3$c.coef) 
  trA.info <- t((wEy *wEyhat)/weights)%*%rep(1,N)
  trA.est <- mean(trA.info)

# find the GCV function
  GCV=  (sum(weights*(residuals)^2)/N )  /( 1- trA.est/np)^2

  
  obj<-list(x=x,y=y,weights=weights, Z=Z,
              d.coef=out1$d.coef, c.coef=out1$c.coef,
              fitted.values=fitted.values, residuals= residuals,
              alpha=alpha, beta=beta, MRinfo=MRinfo,
              overlap=overlap,
              GCV=GCV, lnProfileLike= out2$lnProfileLike,
              rho.MLE=out2$rho.MLE, shat.MLE= out2$shat.MLE,lambda=lambda,
              lnDetCov= out2$lnDetCov, quad.form= out2$quad.form,
              Mc=Mc,
              trA.info=trA.info, trA.est=trA.est,
              eff.df= trA.est, np=np, lambda.fixed=lambda,
              nonzero.entries=nzero,m=2,nt=nt, nZ=nZ,ind.drift=ind.drift,
              call=match.call())
  class(obj)<- "LatticeKrig"
  return(obj)
}

LatticeKrig.coef<- function(Mc,wPHI,wT.matrix,wy,lambda){
  if(length(lambda)>1){
    stop("lambda must be a scalar")}
  A <- forwardsolve(Mc, transpose = TRUE, t(wPHI)%*%wT.matrix, upper.tri = TRUE)
  A <- backsolve(Mc, A)
  A<- t(wT.matrix)%*%(wT.matrix - wPHI%*%A)/lambda
#   A is  (T^t M^{-1} T) 
  b <- forwardsolve(Mc, transpose = TRUE, t(wPHI)%*%wy, upper.tri = TRUE)
  b <- backsolve(Mc, b)  
  b<- t(wT.matrix)%*%(wy- wPHI%*%b)/lambda
# b is   (T^t M^{-1} y)
# Save the intermediate matrix   (T^t M^{-1} T) ^{-1}
# this the GLS covariance matrix of estimated coefficients
# should be small -- the default is 3X3  
  Omega<- solve( A)
# GLS estimates   
  d.coef<- Omega%*% b
# coefficients of basis functions.   
  c.coef <- forwardsolve(Mc, transpose = TRUE, t(wPHI)%*%(wy-wT.matrix%*%d.coef),
                         upper.tri = TRUE)
  c.coef <- backsolve(Mc, c.coef)

  return(list( c.coef=c.coef, d.coef= d.coef))  
  }

LatticeKrig.lnPlike<- function(Mc,Q,y,lambda,residuals,weights){
  N<- length(y)
  Ntotal<- dim(Q)[1]
  c.mKrig<- weights*residuals/lambda
# find log determinant of reg matrix temp for use in the log likeihood
  lnDetReg <- 2 * sum(log(diag(Mc)))
  lnDetQ<-  2* sum( log( diag( chol(Q))))
# now apply a miraculous determinant identity (Sylvester''s theorem)
#  det( I + UV) = det( I + VU)    providing UV is square
# or more generally   
#  det( A + UV) = det(A) det( I + V A^{-1}U)
#  or as we use it
#  ln(det( I + V A^{-1}U)) = ln( det(  A + UV)) - ln( det(A))
#  
#  with A = lambda* t(H)%*%H , U= t(PHI) V= PHI   
#  M =  lambda*(I + V A^{-1}U)
#
#   lnDet (M) = N* ln(lambda) + ln(Det(I + V A^{-1}U))
#             = N* ln(lambda) + ln(Det( temp)) - ln Det( lambda*Q)
#             = N* ln(lambda) + ln(Det( temp))  - ln Det( Q) - Ntotal* log( lambda)
# Here we  canceled out lambda terms 
# 
  lnDetCov<- lnDetReg - lnDetQ + (N - Ntotal)* log(lambda)
# finding quadratic form
# this uses a shortcut formula for the quadratic form in terms of
# the residuals 
  quad.form<-   sum(y* c.mKrig )
# MLE estimate of rho and sigma
# these are derived by assuming Y is  MN(  Td, rho*M )  
  rho.MLE <- quad.form/N
  shat.MLE <- sigma.MLE <- sqrt(lambda * rho.MLE)
# the  log profile likehood with  rhohat  and  dhat substituted
# leaving a profile for just lambda.
# note that this is _not_  -2*loglike just the log and
# includes the constants
  lnProfileLike <- (-N/2 - log(2*pi)*(N/2)
                      - (N/2)*log(rho.MLE) - (1/2) * lnDetCov)
  return( list(lnProfileLike=lnProfileLike,rho.MLE=rho.MLE, shat.MLE=shat.MLE,
               quad.form=quad.form, lnDetCov=lnDetCov) )
 
}


surface.LatticeKrig<- function(obj,...){
  surface.Krig( obj,...)}

predict.LatticeKrig<- function( object, xnew=NULL,Z=NULL,drop.Z=FALSE,...){
  if( is.null(xnew)){
    xnew<- object$x}
  if( is.null(Z)& object$nZ>0){
    Z<- object$Z} 
  NG<- nrow( xnew)
  if( drop.Z|object$nZ==0){
     temp1<-cbind( rep(1,NG), xnew)%*% object$d.coef[object$ind.drift,]}
   else{
     temp1<- cbind( rep(1,NG), xnew,Z)%*% object$d.coef}
    PHIg<-   MRbasis( xnew,object$MRinfo)
    temp2<- PHIg%*%object$c.coef
return( temp1 + temp2)
}

print.LatticeKrig <- function(x, digits=4, ...){
   MRinfo<- x$MRinfo
   if (is.matrix(x$residuals)) {
        n <- nrow(x$residuals)
        NData <- ncol(x$residuals)
    }
    else {
        n <- length(x$residuals)
        NData <- 1
    }
    c1 <- "Number of Observations:"
    c2 <- n
    if (NData > 1) {
        c1 <- c(c1, "Number of data sets fit:")
        c2 <- c(c2, NData)
    }
    c1 <- c(c1, "Degree of polynomial null space ( base model):")
    c2 <- c(c2, x$m - 1)
    c1 <- c(c1, "Number of parameters in the null space")
    c2 <- c(c2, x$nt)
    if( x$nZ>0){
     c1 <- c(c1, "Number of covariates")
     c2 <- c(c2, x$nZ)}
    if (!is.na(x$eff.df)) {
        c1 <- c(c1, " Eff. degrees of freedom")
        c2 <- c(c2, signif(x$eff.df, digits))
        if (length(x$trA.info) < x$np) {
            c1 <- c(c1, "   Standard Error of estimate: ")
            c2 <- c(c2, signif(sd(x$trA.info)/sqrt(length(x$trA.info)), 
                digits))
        }
    }
    c1 <- c(c1, "Smoothing parameter")
    c2 <- c(c2, signif(x$lambda.fixed, digits))
    if (NData == 1) {
        c1 <- c(c1, "MLE sigma ")
        c2 <- c(c2, signif(x$shat.MLE, digits))
        c1 <- c(c1, "MLE rho")
        c2 <- c(c2, signif(x$rho.MLE, digits))
    }
    c1 <- c(c1, "Nonzero entries in covariance")
    c2 <- c(c2, x$nonzero.entries)
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    cat(" ", fill = TRUE)
    cat("Covariance Model: Wendland/Lattice", fill = TRUE)
    if( MRinfo$eflag){
      cat("Exponential tapered basis comprises finest level", fill=TRUE)}
    cat(MRinfo$nlevel, "level(s)", MRinfo$Ntotal, " basis functions", fill=TRUE)
    for( k in 1: MRinfo$nlevel){ 
      cat("Lattice level", k,"is ", MRinfo$N1[k],"X",MRinfo$N2[k], fill=TRUE)}
    cat( "total number of points: ", x$np,  "overlap of ", x$overlap, fill=TRUE)
    cat("Value(s) for weighting (alpha): ", x$alpha,fill = TRUE)
    cat("Value(s) for lattice dependence (beta): ", x$beta,fill = TRUE)
    invisible(x)
}                      
