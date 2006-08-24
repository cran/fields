"predict.Krig" <-
function (object, x = NULL, lambda = NA, df = NA, model = NA, 
    eval.correlation.model = TRUE, y = NULL, yM= NULL, verbose = FALSE, ...) 
{

# y is full data yM are the data collapsed to replicate means

# if new data is not passed then copy from the object

    if (is.null(y)&is.null( yM)) {
        temp.c <- object$c
        temp.d <- object$d
    }

# default is to predict at data x's
    if (is.null(x)) {
        x <- object$x
    }
    x <- as.matrix(x)
 
    if (verbose) 
        print(x)

# transformations of x values used in Krig

    xc <- object$transform$x.center
    xs <- object$transform$x.scale

    x <- scale(x, xc, xs)

# NOTE knots are already scaceld in Krig object and are used 
# in transformed scale. 

#    knots <- scale( object$knots, xc, xs)

#
# figure out if the coefficients for the surface needto be recomputed. 
    find.coef <- (!is.null(y) | !is.null(yM) | !is.na(lambda) 
                          | !is.na(df)|!is.na(model[1]))

    if( verbose){ cat("find.coef", find.coef, fill=TRUE)}

#   convert effective degrees of freedom to equivalent lambda
    if (!is.na(df)) {
        lambda <- Krig.df.to.lambda(df, object$matrices$D)
    }

    if (!is.na(model)) {
        lambda <- model[1]
    }

    if (is.na(lambda)) 
        lambda <- object$lambda

# 
# if the coefficients need to be recomputed  do it. 
    if (find.coef) {

        if( verbose){ cat("new coefs found", fill=TRUE)}

        object3 <- Krig.coef( object, lambda = lambda, y=y, yM=yM)
        temp.d <- object3$d
        temp.c <- object3$c
    
    }

        if (verbose) {
            cat(" d coefs",fill=TRUE)
            print(temp.d)
            cat("c coefs", fill=TRUE)
            print(temp.c)
        }

# 
# The covariance function is 
# evaluated by using it name, do.call function and any 
# additional arguments. 
#
#
# this is the fixed part of predictor 
#
     temp <- c( object$make.tmatrix(x, object$m) %*% temp.d )
#
# now find sum of covariance functions times coefficients
# Note that the multiplication of the cross covariance matrix
# by the coefficients is done implicitly in the covariance function
# 
     temp<- temp + 
      do.call(
        object$cov.function.name, 
       c( object$args, list(x1 = x, x2 = object$knots, C = temp.c)))

# coerce to vector      
       temp<- c( temp)
    
#
# transform back into raw scale if this is a correlation model.
# if y's are in the scale of correlations 
# if so scale by sd and add back in mean 

    correlation.model <- ( object$correlation.model & eval.correlation.model)

    if (correlation.model) {

      if( !is.na( object$sd.obj[1])){
       temp<- temp*predict( object$sd.obj, x)}

      if( !is.na( object$mean.obj[1])){
       temp<- temp + predict( object$mean.obj, x)}

    }

        return(temp)
}

