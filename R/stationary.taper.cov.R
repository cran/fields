"stationary.taper.cov" <-
 function (x1, x2, Covariance="Exponential", 
           Distance="rdist",Taper="Wendland", 
           Dist.args=NULL, Taper.args=NULL, 
           theta=1.0,C=NA, marginal=FALSE,
           max.points=NULL, mean.neighbor=50,
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
        else{
                Taper.args<- c( Taper.args, 
                               list(dimension=ncol(x1))) }
     }  

#
# Add in general defaults for taper arguments if not Wendland
#  theta = 1.0 is the default range for the taper. 

     if( is.null(Taper.args)){
               Taper.args<- list( theta=1.0) } 

# default for max.points work array -- see help(rdist) for details)

    if( is.null( max.points) ) { max.points<- nrow( x1)*mean.neighbor }


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
# find nearest neighbor distances based on taper threshhold

    fields.rdist.near( x1,x2,delta= taper.range, max.points= max.points,
                          mean.neighbor=mean.neighbor)-> sM
    
# ra components are the pairwise distances
# apply covariance and taper to these.     
     sM$ra<- (do.call(Covariance,  c(list(d=sM$ra),Cov.args) )*
             do.call(Taper,  c(list(d=sM$ra),Taper.args) ) )

# at this point sM is in index format with the component ind 
# giving row and column indices and ra being the nonzero values.


# decide whether to return sM in spam sparse form or as a full matrix
#
# SPAM format is column index (ja) and the position in the array for the beginning
# of each row (ia) see spind2spam for details.
# 

 if( spam.format){ 
     return(spind2spam(sM)) }
 else{
     return( spind2full(sM)) }
}


# or multiply cross covariance by C
 if( !is.na( C[1]) ){
# 
# as coded below this is not particularly efficient of memory 
# 
# see comments in previous if block for this code
   fields.rdist.near( x1,x2,delta= taper.range, max.points= max.points)-> sM
   sM$ra<- (do.call(Covariance,  c(list(d=sM$ra),Cov.args) )*
             do.call(Taper,  c(list(d=sM$ra),Taper.args) ) )

   if( spam.format) {
       return(  spind2spam(sM) %*% C )}
   else{
       return( spind2full(sM) %*% C)  }
 }


# or find marginal variance and return  a vector. 
  if( marginal){
    sigma2<- do.call(Covariance, c(list(d=0),Cov.args) ) *  do.call(Taper,  c(list(d=0),Taper.args) )
    return( rep( sigma2, nrow( x1)))}  

# should not get here!
}

