# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"stationary.taper.cov" <-
 function (x1, x2, Covariance="Exponential",
           Taper="Wendland",
           Dist.args=NULL, Taper.args=NULL,
           theta=1.0,C=NA, marginal=FALSE,
           spam.format=TRUE, ... )
{


# get covariance function arguments from call

 Cov.args<- list(...)

# coerce x1 and x2 to matrices

    if( is.data.frame( x1)) x1<- as.matrix(x1)

    if (!is.matrix(x1))
        x1 <- matrix(c(x1), ncol=1)

    if (missing(x2))
        x2 <- x1

    if( is.data.frame( x2)) x2<- as.matrix(x1)

    if (!is.matrix(x2))
        x2 <- matrix(c(x2), ncol=1)

# Default taper arguments that are particular to the Wendland.
     if( Taper=="Wendland"){
        if(is.null(Taper.args)) {
                Taper.args<- list(
                    theta=1.0, k=2, dimension=ncol(x1)) }
#rf        else{
#rf                Taper.args<- c( Taper.args,
#rf                               list(dimension=ncol(x1))) }
     } 

#
# Add in general defaults for taper arguments if not Wendland
#  theta = 1.0 is the default range for the taper.

     if( is.null(Taper.args)){
               Taper.args<- list( theta=1.0) }


#########################################
# considering the scale parameters

    if (length(theta) == 1)
        theta <- rep(theta, ncol(x1))
# handle special case of 1-d
    if( ncol(x1)==1) { theta<- matrix( c(theta),1,1)}
# handle special case of just diagonal elements of  theta
    if (is.vector(theta))
        theta <- diag(theta)
# following now treats theta as a full matrix for scaling and rotation.
#
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    x1 <- x1 %*% t(solve(theta))
    x2 <- x2 %*% t(solve(theta))
#
# locations are now scaled and rotated correctly

# copy taper range

      taper.range<- Taper.args$theta

#NOTE tapering is applied to the _scaled_ locations.

# now apply covariance function to pairwise distance matrix, or multiply
# by C vector or just find marginal variance


 if(  !marginal & is.na(C[1]) ){

# this if block finds the cross covariance matrix
# find nearest neighbor distances based on taper threshhold.
# This is hardwired to 'nearest.dist' function from spam.    

   sM <- do.call('nearest.dist',c(list(x1, x2, delta = taper.range,
                                       upper = NULL, diag=TRUE), Dist.args))
# sM@entries contains are the pairwise distances
# apply covariance and taper to these.
   sM@entries <- do.call(Covariance, c(list(d = sM@entries), Cov.args)) *
                 do.call(Taper, c(list(d = sM@entries), Taper.args))


# decide whether to return sM in spam sparse form or as a full matrix

    if (spam.format) {
      return(sM)
    }
    else {
      return(as.matrix(sM))
    }

 }  
# or multiply cross covariance by C
 if( !is.na( C[1]) ){
#
# see comments in previous if block for this code
   sM <- do.call('nearest.dist',c(list(x1, x2, delta = taper.range,
                                       upper = NULL, diag=TRUE), Dist.args))
   sM@entries <- do.call(Covariance, c(list(d = sM@entries), Cov.args)) *
                 do.call(Taper, c(list(d = sM@entries), Taper.args))
   
   return(sM %*% C)
 }


# or find marginal variance and return  a vector.
  if( marginal){
    sigma2 <- do.call(Covariance, c(list(d=0),Cov.args) ) *
              do.call(Taper,  c(list(d=0),Taper.args) )
    return( rep( sigma2, nrow( x1)))}

# should not get here!
 }

