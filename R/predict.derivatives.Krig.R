# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"predict.derivative"<- function( object,...){
   UseMethod("predict.derivative")}  

"predict.derivative.Krig" <- function(object, x = NULL, 
    verbose = FALSE, ...) {
    # this is a lean evaluation of the derivatives of the
    # random component of the model.
    # several checks to make sure this being applied to
    # simple Krig models where it makes sense
    if (object$correlation.model) {
        stop("Can not handle correlation model with derivative evaluation")
    }
    if (object$null.function.name != "Krig.null.function") {
        stop("null space may not be a low order polynomial")
    }
    # default is to predict at data x's
    if (is.null(x)) {
        x <- object$x
    }
    else {
        x <- as.matrix(x)
    }
    # transformations of x values used in Krig
    xc <- object$transform$x.center
    xs <- object$transform$x.scale
    x <- scale(x, xc, xs)
    # NOTE knots are already scaled in Krig object and are used
    # in transformed scale.
    #  i.e.   knots <- scale( object$knots, xc, xs)
    temp.d <- object$d
    temp.c <- object$c
    if (verbose) {
        cat(" d coefs", fill = TRUE)
        print(temp.d)
        cat("c coefs", fill = TRUE)
        print(temp.c)
    }
    #
    # this is the polynomial fixed part of predictor
    #
    temp1 <- fields.derivative.poly(x, m = object$m, object$d)
    # add in spatial piece
    # The covariance function is
    # evaluated by using it name, do.call function and any
    # additional arguments.  Note use of derivative and 'C' arguments
    # to do multiplication of partials of covariance times the C
    # vector. If C is a matrix of coefficients a error is produced.
    temp2 <- do.call(object$cov.function.name, c(object$args, 
        list(x1 = x, x2 = object$knots, derivative = 1, C = temp.c)))
    # returned value is the matrix of partials of polynomial plus  partials of spatial # part aso add in chain rule scale factor  because
    # functional form for the surface uses the coordinates xscaled =  (x- xc)/xs
      return(t(t(temp1 + temp2)/xs))
}
