Krig.check.xY<- function( x,Y, weights, na.rm,verbose=FALSE){

#
# check for missing values in Y or X.
#
# save logical indicating where there are NA's  
# and check for NA's    
#
    ind<- is.na( Y)
    if (any(ind) & !na.rm) {
        stop("Need to remove missing values or use: na.rm=TRUE in the call")
    }

#
# coerce x to be a matrix
    x <- as.matrix(x)
#
# coerce Y to be a vector
#
   Y <- c(Y)


#
#default weights ( reciprocal variance of errors). 
#
    if (is.null(weights)) 
       weights <- rep(1, length(Y))
#
# check that dimensions agree
#
    if(length( Y) != nrow( x)){
       stop(" length of y and number of rows of x differ" )}
    if(length( Y) != length(weights)){
       stop(" length of y and weights differ" )}


#
# if NAs can be removed then remove them and warn the user
    if (na.rm) {
        ind <- is.na(Y)
        if (any(ind)) {
            Y <- Y[!ind]
            x <- x[!ind, ]
            weights <- weights[!ind]
#            warning("NA's have been removed from Y ")
        }
    }
#
# check for NA's in x matrix -- there should not be any to proceed!

    if( any( c( is.na(x)))){ stop(" NA's in x matrix")}

#
# verbose block 
    if (verbose) {
        cat("Y:",fill=TRUE)
        print(Y)
        cat("x:",fill=TRUE)
        print(x)
        cat("weights:", fill=TRUE)
        cat(weights, fill=TRUE)
    }

#
# save x, weights  and Y w/o NAs
    N <- length(Y)
   return( list(N= N, y=Y, x=x, weights=weights, NA.ind= ind))

}

