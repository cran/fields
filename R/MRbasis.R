# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

MRbasis<- function(x1, MRinfo, spam.format = TRUE, verbose = FALSE){
# order of Wendland is hardwired
  Korder<- 2
  grid.info<- MRinfo$grid.info
  nlevel<- MRinfo$nlevel
  eflag<- MRinfo$eflag
# INFLATE and overlap are hardwired in LatticeKrig  
  INFLATE<- MRinfo$INFLATE
  overlap<-MRinfo$overlap
#  
# accumulate matrix column by column in PHI  
  PHI<-NULL
  for( j in 1:nlevel){
# loop over levels the last level might be special ...      
      delta<- MRinfo$delta[j]
      grid.list<- list( x= seq( grid.info$xmin, grid.info$xmax,delta),
                        y= seq( grid.info$ymin, grid.info$ymax,delta))
      centers<- make.surface.grid( grid.list)
      if(verbose){
        print(dim(centers))}
#  set the range of basis functions, they are assumed to be zero outside
#  the radius basis.delta
      basis.delta<- delta*overlap 
# finest level can be different: if  eflag==TRUE uase a more direct Kriging basis
# such as an exponential.
      if( !(( j==nlevel)&(eflag)) ){
      PHI<-cbind( PHI,
                  wendland.cov( x1,centers, theta= basis.delta, k=Korder,
                                             spam.format=spam.format))}
      else{
#        cat("exponent basis functions","level", j, fill=TRUE)
# columns index basis functions so the different levels are accumulated
# column by column       
      PHI<- cbind(PHI,
                   stationary.taper.cov( x1,centers, theta= basis.delta,
                     spam.format=spam.format,                                                                Taper.args=list( theta= basis.delta*INFLATE,                                                               k=Korder, dimension=2)))}   
  }
# attach  MRinfo list to the matrix to help identify how the basis functions
# are organized. 
  attr(PHI, which="info")<-MRinfo
  return(PHI)
}


  MRbasis.setup<- function(grid.info, nlevel=1, eflag=FALSE,overlap=2.5, INFLATE=1.5){
# determines the multiresolution basis function indices.
#     
    delta<- grid.info$delta  
    grid.list<- list( x= seq( grid.info$xmin, grid.info$xmax,delta),
                    y= seq( grid.info$ymin, grid.info$ymax,delta))
  
    delta.save<- N1<-N2<- rep(NA, nlevel)
    N1[1]<- length( grid.list$x)
    N2[1]<- length( grid.list$y)
    delta.save[1]<- delta
 
    if( nlevel >1){
      for( j in 2:nlevel){     
        delta<- delta/2
        delta.save[j]<- delta
        grid.list<- list( x= seq( grid.info$xmin, grid.info$xmax,delta),
                    y= seq( grid.info$ymin, grid.info$ymax,delta))
        N1[j]<- length( grid.list$x)
        N2[j]<- length( grid.list$y)      
      }    
    }
    offset<- c( 0, cumsum(N1*N2))
    list( N1=N1, N2=N2,nlevel=nlevel,delta= delta.save,Ntotal= sum( N1*N2),
                 offset=offset,grid.info= grid.info, eflag=eflag,
                 overlap=overlap,INFLATE=INFLATE)
}


MRprecision <- function( MRinfo, alpha=1.0, beta=-.2, spam.format = TRUE){
   eflag<- MRinfo$eflag
   N1<- MRinfo$N1
   N2<- MRinfo$N2
   INFLATE<- MRinfo$INFLATE
   grid.info<- MRinfo$grid.info
   NL<- MRinfo$nlevel
   if(NL!= length(N2)){
     stop("N1 and N2 are not equal")}
   offset<- MRinfo$offset
# fill out alpha and beta if just a single value   
   if( length(alpha)==1){   
     alpha<- rep( alpha, NL)}   
   if( length(beta)==1){   
     beta<- rep( beta, NL)}
#   cat("beta", beta, fill=TRUE)
# ind holds non-zero indices and ra holds the values   
   ind<- NULL
   ra<-NULL
# loop over levels the last level might be special ...   
   for( j in 1:NL){
      delta.temp<- MRinfo$delta[j]
# BAU -- business as usual (no exponential fill in)      
      BAU<-  ifelse(j==NL,!eflag,TRUE)
#      print( BAU)
      if(BAU){
#         print(beta[j])
         temp<-lattice.precision(N1[j], N2[j],beta= beta[j], spam.format=TRUE)
         temp<- temp%*%t(temp)}
      else{
# if eflag is set fill in last level with a tapered exponential covariance        
           delta.temp<- MRinfo$delta[j]
           centers<- make.surface.grid(
                       list( x= seq( grid.info$xmin, grid.info$xmax,delta.temp),
                          y= seq( grid.info$ymin, grid.info$ymax,delta.temp)))
           temp<-stationary.taper.cov( centers,centers, theta= delta.temp,
                        spam.format=spam.format,
                        Taper.args=list( theta=delta.temp*INFLATE,k=3, dimension=2))
         }
# convert back to index format to make it easier to accumulate the blocks of Q          
      temp<- spam2spind(temp)
# accumulate the new block in the growing matrix.            
      ind<- rbind( ind, temp$ind+offset[j])
      ra<- c( ra,alpha[j]*temp$ra)   
    }
   da<- c( offset[NL+1], offset[NL+1])
   
   if( spam.format){
     return( spind2spam( list( ind=ind, ra=ra, da=da) ))}
   else{
     return( list( ind=ind, ra=ra, da=da) )}
 }

"lattice.precision" <- function(N1,N2, spam.format = TRUE, beta) {
  NN<- N1*N2
  da<- as.integer( c(NN,NN))
  I<- as.integer(rep(1:N1,N2))
  J<- as.integer(rep((1:N2), rep( N1,N2)))
# contents of sparse matrix
  ra<- rep( 1, N1*N2)
  ra<- c( ra, rep( beta, 4*NN))
#  Note: if beta depends on lattice position
# pass beta as a matrix and use  ra<- as.numeric(c( ra,rep( c(beta), 4)))
  Hi<- rep(1:NN,5)
# if periodic wrap row indices from off of grid to other side
 
  Hj<- c( I    + (J-1)*N1,
         (I-1) + (J-1)*N1,
         (I+1) + (J-1)*N1,
          I    + (J-2)*N1,
          I    +   (J)*N1 )
  good<- c( rep( TRUE,NN),
           (I-1) > 0,
           (I)   < N1,
           (J-2) >=0,
           (J)   <N2 )
# remove cases that are beyond the lattice (edges, corners)        
  Hi<- as.integer(Hi[good])
  Hj<- as.integer(Hj[good])
  ra<- ra[good]
# Convert matrix to spam format  
  sM <- spind2spam( list( ind=cbind(Hi,Hj), ra=ra, da=da) )
# convert to an ordinary matrix (usually for debugging)     
  if (!spam.format) {
          sM<- as.matrix(sM)}
  return(sM)
}


MR.cov<-function( x1,x2=NULL,grid.info, alpha=1.0, nlevel=1,
                          beta=-.2,
                          C=NA, marginal=FALSE, eflag=FALSE, INFLATE=1.5){
  if( is.null(x2)){ x2<-x1}
  MRinfo<- MRbasis.setup( grid.info,  nlevel=nlevel, eflag=eflag, INFLATE=INFLATE)
  PHI1<- MRbasis( x1,MRinfo)
  PHI2<- MRbasis( x2,MRinfo)
  Q<-  MRprecision(MRinfo,alpha=alpha, beta=beta, spam.format=TRUE)
#  print( beta)
#  Note:  Q = (H)%*%t(H)
# find Q^{-1} %*% t(PHI2)  
  chol(Q)-> Mc 
 
   if (is.na(C[1]) & !marginal) {
     A <- forwardsolve(Mc, transpose = TRUE, t(PHI2), upper.tri = TRUE)
     A <- backsolve(Mc, A)
     return( PHI1%*%A)
   }
   if (!is.na(C[1])) {
     A <- forwardsolve(Mc, transpose = TRUE, t(PHI2)%*%C, upper.tri = TRUE)
     A <- backsolve(Mc, A)
        return(PHI1%*%A)
    }
    if (marginal) {
      stop("marginal not implemented")
        return(rep(1, nrow(x1)))
    }


}

MR.sim<- function(x1,grid.info, alpha=1.0, nlevel=1, beta=-.2, overlap=2.5,
                  eflag=FALSE, INFLATE=1.5, M=1){
  MRinfo<- MRbasis.setup( grid.info, nlevel=nlevel, eflag=eflag, INFLATE=1.5)
  PHI1<- MRbasis(x1,MRinfo)
  Q<-  MRprecision(MRinfo,alpha=alpha, beta=beta, spam.format=TRUE)
#  
#  Q is precision matrix of the coefficients -- not of the field
#  last step does the multiplication to go from coefficients to evaluating
#  values at the field
#  Q = t(H)%*%H = inv((Sigma)
#  So   Sigma= Hi%*% t(Hi)  
#  find u= t(Hi) %*% N(0,1)   then cov(u) = t(Hi)%*%Hi
#  Hi is upper triangular
#
# snippet of code to test the algebra ...  
#   x<-seq( 0,1,,20); Sigma<- exp(-rdist( x,x)/2.5); Q<- solve( Sigma)
#   chol(Q)->Mc; H<- Mc ; Hi<- solve(H);
#   test.for.zero( Q, t(H)%*%H); test.for.zero(Sigma, Hi%*%t(Hi))   
#   E<- rnorm(20);  u1<- Hi%*% E ;   u2<-backsolve(Mc,E)
#   test.for.zero(u1,u2)  
# 
  chol(Q)-> Mc
  E<- matrix( rnorm( M*MRinfo$Ntotal), nrow=MRinfo$Ntotal,ncol= M)
  A <- backsolve(Mc,E)
  return( PHI1%*%A)
}
