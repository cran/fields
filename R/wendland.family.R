# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


Wendland2.2 <- function(d, theta=1.0)
{
# Cari's test function  with explicit form  for d=2 k=2
# taper range is 1.0
 d<- d/theta
 if (any(d<0))
   stop("d must be nonnegative")
 return(((1-d)^6*(35*d^2+18*d+3))/3*(d<1))
}

# Tapering function

Wendland<- function( d,theta=1.0, dimension=2,k, derivative=0, phi=1){
    # d = dimension k = order. see wendland.coef for details

#    if( missing( dimension)){ 
#      stop("need to give a dimension argument to Wendland function")}
     
    if( any(d<0) ){
      stop("some distances are negative")}
     coef<- wendland.coef( dimension,k)

#  take derivative of polynomial if derivative > 0
    if( derivative !=0) {
        L<- length(coef)
        coef<- coef[2:L]* (1:(L-1))/theta 
      }
    
   ifelse( d<=theta, phi*fields.evlpoly( d/theta, coef), 0)

  }

# the monster
"wendland.cov" <-
function (x1, x2, theta = rep(1, ncol(x1)), k=2,
          C = NA, marginal=FALSE,
          Dist.args=NULL,
          spam.format=TRUE,
           derivative=0)
{

#
#   if marginal variance is needed
#  this is a quick return

    if( marginal){
    return( rep( 1, nrow(x1)) )}

#  the rest of the possiblities require some computing

     if (!is.matrix(x1))
        x1 <- as.matrix(x1)
    if (missing(x2))
        x2 <- x1
    if (!is.matrix(x2))
        x2 <- as.matrix(x2)
    if (length(theta) == 1)
        theta <- rep(theta, ncol(x1))
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
# scale the coordinates by theta
# a more general scaling by a matrix is done in stationary.cov
    x1 <- scale(x1 , center=FALSE, scale=theta)
    x2 <- scale(x2 , center=FALSE, scale=theta)

# once scaling is done taper is applied with default range of 1.0

# find polynomial coeffients that define
# wendland on [0,1]
#  d dimension and  k is the order



#  find sparse matrix of Euclidean distances
#  ||x1-x2||**2 (or any other distance that may be specified by
#  Dist.args    

   sM <- do.call('nearest.dist',c(list(x1, x2, delta = 1, upper = NULL,
                     diag=TRUE), Dist.args))

#    diag(sM) <- 0    # required as passing two matrices

# equivalent to sM<- rdist( x1,x2) but sM is in sparse index format


#
# there are two possible actions listed below:

# find cross covariance matrix
# return either in sparse or matrix format


    if (is.na(C[1]) ) {

       sM@entries<- Wendland( sM@entries, theta=1.0,k=k, dimension=d)

       if(!spam.format){
                 return( as.matrix(sM) )}
            else{
                 return(  sM )}
    }
     else{
#
# multiply cross covariance matrix by C
#  note multiply happens in spam matrix format
#


        if( derivative ==0){
           sM@entries <- Wendland( sM@entries,dimension=d, theta=1.0,k=k)
          return( sM%*%C) }

#        otherwise evaluate partial derivatives with respect to x1
#        Big mess of code and an explicit for loop!

        else {

          L<- length( coef)
#         loop over dimensions and accumulate partial derivative matrix.

          tempD<- sM@entries           # save nonzero  distances

          tempW<- Wendland( tempD, theta=1.0,k=k,dimension=d,
                                            derivative=derivative)

# loop over dimensions and knock out each partial accumalte these in
# in  temp

          temp<- NULL
# Create rowindices vector          
          sMrowindices <- rep(1:n1, diff(sM@rowpointers))
          for( kd in 1:d){
#
#            Be careful if the distance (tempD) is close to zero.
#            Also adjust for the fact that the x1 and x2 are scaled


             sM@entries <- ifelse( tempD== 0, 0,
              ( x1[sMrowindices,kd]- x2[sM@colindices,kd] )/ (theta[kd]*tempD)*tempW )

# accumlate the new partial
             temp<- cbind( temp, sM%*%C )  }

        return( temp)}

   } # end of  C evaluation block


# should not get here!
  }



`wendland.coef` <-
function(d,k){

# Calculate coefficients according to Thm 1.3 of Wendland (1998)
 l <- floor(d/2) + k + 1
 coef <- matrix(NA, nrow=l+2*k+1, ncol=k+1)
 coef[1:(l+1),1] <- (-1)^(0:l) * choose(l, 0:l)

 if(k>0)
   for(s in 0:(k-1)){
     i <- seq(0,l+2*s)
     coef[1,s+2] <- sum(coef[i+1,s+1]/(i+2))
     coef[2,s+2] <- 0
     for(j in 2:(l+2*s+2)){
       coef[j+1,s+2] <- -coef[j-1,s+1]/j
     }
   }
 final.coef <- coef[,k+1]/coef[1,k+1]

return( final.coef)
}
