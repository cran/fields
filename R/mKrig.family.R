# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

mKrig<- function (x,y,
                  weights = rep( 1, nrow( x)), lambda=0,
                  cov.function= "stationary.cov", m=2,
                  chol.args = NULL, cov.args=NULL,...)
{

# grab other arguments for covariance function
    cov.args<- c( cov.args, list( ...))
   

# see comments in Krig.engine.fixed for algorithmic commentary
#    
#
# default values for Cholesky decomposition, these are important
# for sparse matrix decompositions 

    
# check for duplicate x's.
# stop if there are any
        if( any(duplicated( cat.matrix(x))) ) 
                  stop("locations are not unique see help(mKrig) ")

# create fixed part of model as m-1 order polynomial 
         Tmatrix <- fields.mkpoly(x, m)
        
# set some dimensions
         np <- nrow( x)
         nt <- ncol( Tmatrix)

# covariance matrix at observation locations
# NOTE: if cov.function is a sparse constuct then tempM will be sparse.
# see e.g. wendland.cov

         tempM<-  do.call(
             cov.function, c(cov.args, list(x1 = x, x2 = x))  )
#
# decide how to handle the pivoting.
# one wants to do pivoting if the matrix is sparse. 
# if tempM is not a matrix assume that it is in sparse format. 
#
    sparse.flag<- !is.matrix( tempM)

#
# set arguments that are passed to cholesky
#
    if( is.null( chol.args)) {
             chol.args<- list( pivot= sparse.flag)}
    else{
         chol.args<- chol.args}

# record sparsity of tempM 

        nzero <- ifelse( sparse.flag, length( tempM@entries), np^2) 

# add diagonal matrix that is the observation error Variance 
# NOTE: diag must be a overloaded function to handle  sparse format.

      if( lambda != 0) {
        diag(tempM) <-   (lambda/weights) + diag(tempM)
    
      }

# cholesky decoposition of tempM
# do.call used to supply other arguments to the function
# especially for sparse applications.
# If chol.args is NULL then this is the same as
#              Mc<-chol(tempM)
    
     Mc<-  do.call("chol",c( list( x=tempM), chol.args) ) 

# Efficent way to multply inverse of Mc times the Tmatrix
    
 VT<- forwardsolve( Mc, x=Tmatrix, transpose=TRUE, upper.tri=TRUE)
 qr( VT)-> qr.VT

# start linear algebra to find solution 
# Note that all these expressions make sense if y is a matrix
# of several data sets and one is solving for the coefficients 
# of all of these at once. In this case d.coef and c.coef are matrices

#
# now do generalized least squares for d 

 d.coef<- qr.coef( qr.VT,
                  forwardsolve( Mc, transpose=TRUE, y, upper.tri=TRUE) )

# and then find c. 
# find the coefficents for the spatial part.    
 
 c.coef<- forwardsolve( Mc, transpose=TRUE,
                        y - Tmatrix%*% d.coef, upper.tri=TRUE)
 c.coef<- backsolve( Mc,c.coef)
    

# return coefficients and   include lambda as a check because 
# results are meaningless for other values of lambda 
# returned list is an "object" of class mKrig (micro Krig)
    
    out<-  list( d=(d.coef), c=(c.coef), 
                nt=nt, 
                np=np,
                lambda.fixed=lambda,
                x= x,
                cov.function.name=cov.function, args=cov.args,m=m,
                chol.args=chol.args,
                call= match.call(), 
                nonzero.entries= nzero)
#       
   out$fitted.values<- predict.mKrig( out)
   out$residuals<-  y - out$fitted.values
   

class(out)<- "mKrig"    

return(out)

  }

print.mKrig<- function(x,...){

    c1 <- "Number of Observations:"
    c2 <- length(x$residuals)
  
        c1 <- c(c1, "Degree of polynomial null space ( base model):")
        c2 <- c(c2, x$m - 1)

    c1 <- c(c1, "Number of parameters in the null space")
    c2 <- c(c2, x$nt)

    c1 <- c(c1, "Smoothing parameter")
    c2 <- c(c2, x$lambda.fixed)

    c1 <- c(c1, "Nonzero entries in covariance")
    c2 <- c(c2, x$nonzero.entries)
     
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    cat("Covariance Model:", x$cov.function, fill = TRUE)


    if (x$cov.function == "stationary.cov") {
        cat("  Covariance function is ", x$args$Covariance, fill = TRUE)
    }
    if (!is.null(x$args)) {
        cat("  Names of non-default covariance arguments: ", 
            fill = TRUE)
        cat("      ", paste(as.character(names(x$args)), collapse = ", "), 
            fill = TRUE)
    }

    invisible(x)
  }    




predict.mKrig<- function(object, xnew=NULL, derivative=0,...){


# the main reason to pass new args to the covariance is to increase
# the temp space size for sparse multiplications
# other optional arguments from mKrig are passed along in the
# list object$args

   cov.args<- list(...)

# predict at observation locations by default
  
  if( is.null( xnew)){
        xnew<- object$x }

# fixed part of the model this a polynomial of degree m-1    
# Tmatrix <- fields.mkpoly(xnew, m=object$m)
# 

  if( derivative==0){
  temp1 <- fields.mkpoly(xnew, m=object$m) %*% object$d}
  else{
       temp1<- fields.derivative.poly( xnew,m=object$m, object$d)}

# add nonparametric part. Covariance basis functions
# times coefficients.
# syntax is the name of the function and then a list with
# all the arguments. This allows for different covariance functions
# that have been passed as their name.


  if( derivative == 0){
# argument list are the parameters and other options from mKrig
#  locations and coefficients,
  temp2 <-  do.call(object$cov.function.name, 
             c(object$args,
             list(x1 = xnew, x2 = object$x, C = object$c), cov.args ) )}
  else{
  temp2 <-  do.call(object$cov.function.name, 
             c(
                object$args,
                list(x1 = xnew, x2 = object$x, C = object$c,
                           derivative=derivative), 
                cov.args )
            )}



# add two parts together and coerce to vector 
 return( (temp1 + temp2) )

}


