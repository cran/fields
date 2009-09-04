# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
mKrig <- function(x, y, weights = rep(1, nrow(x)), 
    lambda = 0, cov.function = "stationary.cov", m = 2, chol.args = NULL, 
    cov.args = NULL, find.trA = TRUE, NtrA = 20, iseed = 123, 
    ...) {
    # grab other arguments for covariance function
    cov.args <- c(cov.args, list(...))
    # see comments in Krig.engine.fixed for algorithmic commentary
    #
    #
    # default values for Cholesky decomposition, these are important
    # for sparse matrix decompositions
    # check for duplicate x's.
    # stop if there are any
    if (any(duplicated(cat.matrix(x)))) 
        stop("locations are not unique see help(mKrig) ")
    if (any(is.na(y))) 
        stop("Missing values in y should be removed")
    # create fixed part of model as m-1 order polynomial
    Tmatrix <- fields.mkpoly(x, m)
    # set some dimensions
    np <- nrow(x)
    nt <- ncol(Tmatrix)
    # as a place holder for reduced rank Kriging, distinguish between
    # observations locations and locations to evaluate covariance.
    # (this is will also predict.mKrig to handle a Krig object)
    knots <- x
    # covariance matrix at observation locations
    # NOTE: if cov.function is a sparse constuct then tempM will be sparse.
    # see e.g. wendland.cov
    tempM <- do.call(cov.function, c(cov.args, list(x1 = x, x2 = x)))
    #
    # decide how to handle the pivoting.
    # one wants to do pivoting if the matrix is sparse.
    # if tempM is not a matrix assume that it is in sparse format.
    #
    sparse.flag <- !is.matrix(tempM)
    #
    # set arguments that are passed to cholesky
    #
    if (is.null(chol.args)) {
        chol.args <- list(pivot = sparse.flag)
    }
    else {
        chol.args <- chol.args
    }
    # record sparsity of tempM
    nzero <- ifelse(sparse.flag, length(tempM@entries), np^2)
    # add diagonal matrix that is the observation error Variance
    # NOTE: diag must be a overloaded function to handle  sparse format.
    if (lambda != 0) {
        diag(tempM) <- (lambda/weights) + diag(tempM)
    }
    # cholesky decoposition of tempM
    # do.call used to supply other arguments to the function
    # especially for sparse applications.
    # If chol.args is NULL then this is the same as
    #              Mc<-chol(tempM)
    Mc <- do.call("chol", c(list(x = tempM), chol.args))
    lnDetCov <- 2 * sum(log(diag(Mc)))
    # Efficent way to multply inverse of Mc times the Tmatrix
    VT <- forwardsolve(Mc, x = Tmatrix, transpose = TRUE, upper.tri = TRUE)
    qr.VT <- qr(VT)
    # start linear algebra to find solution
    # Note that all these expressions make sense if y is a matrix
    # of several data sets and one is solving for the coefficients
    # of all of these at once. In this case d.coef and c.coef are matrices
    #
    # now do generalized least squares for d
    d.coef <- qr.coef(qr.VT, forwardsolve(Mc, transpose = TRUE, 
        y, upper.tri = TRUE))
    # and then find c.
    # find the coefficents for the spatial part.
    c.coef <- forwardsolve(Mc, transpose = TRUE, y - Tmatrix %*% 
        d.coef, upper.tri = TRUE)
    # save intermediate result this is   t(y- T d.coef)( M^{-1}) ( y- T d.coef)
    quad.form <- sum(c.coef^2)
    # find c coefficients
    c.coef <- backsolve(Mc, c.coef)
    # GLS covariance matrix for fixed part.
    Rinv <- solve(qr.R(qr.VT))
    Omega <- Rinv %*% t(Rinv)
    # MLE estimate of rho and sigma
    rhohat <- sum(c.coef * y)/(np - nt)
    rho.MLE <- quad.form/np
    shat.MLE <- sigma.MLE <- sqrt(lambda * rho.MLE)
    # the  log profile likehood with  rhohat  and  dhat substituted
    # leaving a profile for just lambda.
    # note that this is _not_  -2*loglike just the log and
    # includes the constants
    lnProfileLike <- (-np/2 - log(2*pi)*(np/2)
                      - (np/2)*log(rho.MLE) - (1/2) * lnDetCov)
    #
    # return coefficients and   include lambda as a check because
    # results are meaningless for other values of lambda
    # returned list is an 'object' of class mKrig (micro Krig)
    # also save the matrix decopositions so coefficients can be
    # recalculated for new y values.
    out <- list(d = (d.coef), c = (c.coef), nt = nt, np = np, 
        lambda.fixed = lambda, x = x, knots = knots, cov.function.name = cov.function, 
        args = cov.args, m = m, chol.args = chol.args, call = match.call(), 
        nonzero.entries = nzero, shat.MLE = shat.MLE, rho.MLE = rho.MLE, 
        sigma.MLE, rhohat = rhohat, lnProfileLike = lnProfileLike, 
        lnDetCov = lnDetCov,quad.form = quad.form, Omega = Omega, qr.VT = qr.VT, Mc = Mc, 
        Tmatrix = Tmatrix)
    #
    # find the residuals directly from solution 
    # to avoid a call to predict 
    out$residuals<- -lambda*c.coef
    out$fitted.values <- y - out$residuals
    # estimate effective degrees of freedom using Monte Carlo trace method.
    if (find.trA) {
        out2 <- mKrig.trace(out, iseed, NtrA)
        out$eff.df <- out2$eff.df
        out$trA.info <- out2$trA.info
        out$GCV <- (sum(out$residuals^2)/np)/(1 - out2$eff.df/np)^2
        if (NtrA < np) {
            out$GCV.info <- (sum(out$residuals^2)/np)/(1 - out2$trA.info/np)^2
        }
        else {
            out$GCV.info <- NA
        }
    }
    else {
        out$eff.df <- NA
        out$trA.info <- NA
        out$GCV <- NA
    }
    class(out) <- "mKrig"
    return(out)
}
mKrig.trace <- function(object, iseed, NtrA) {
    set.seed(iseed)
    # if more tests that number of data points just
    # find A exactly by np predicts.
    if (NtrA >= object$np) {
        Ey <- diag(1, object$np)
        NtrA <- object$np
        trA.info <- diag(predict.mKrig(object, ynew = Ey))
        trA.est <- sum(trA.info)
    }
    else {
        # if fewer tests then use random trace method
        # find fitted.values  for iid N(0,1) 'data' to calculate the
        # the Monte Carlo estimate of tr A(lambda)
        # basically repeat the steps above but take some
        # short cuts because we only need fitted.values
        # create random normal 'data'
        Ey <- matrix(rnorm(object$np * NtrA), nrow = object$np, 
            ncol = NtrA)
        trA.info <- colSums(Ey * (predict.mKrig(object, ynew = Ey)))
        trA.est <- mean(trA.info)
    }
    return(list(trA.info = trA.info, eff.df = trA.est))
}
mKrig.coef <- function(object, y) {
    # given new data y and the matrix decompositions in the
    # mKrig object find coefficients d and c.
    # d are the coefficients for the fixed part
    # in this case hard coded for a low order polynomial
    # c are coefficients for the basis functions derived from the
    # covariance function.
    #
    # see mKrig itself for more comments on the linear algebra
    #
    # Note that all these expressions make sense if y is a matrix
    # of several data sets and one is solving for the coefficients
    # of all of these at once. In this case d.coef and c.coef are matrices
    #
    # generalized least squares for d
    d.coef <- qr.coef(object$qr.VT, forwardsolve(object$Mc, transpose = TRUE, 
        y, upper.tri = TRUE))
    #  residuals from subtracting off fixed part
    #  of model as m-1 order polynomial
    Tmatrix <- fields.mkpoly(object$x, object$m)
    resid <- y - Tmatrix %*% d.coef
    # and now find c.
    c.coef <- forwardsolve(object$Mc, transpose = TRUE, resid, 
        upper.tri = TRUE)
    c.coef <- backsolve(object$Mc, c.coef)
    out <- list(d = (d.coef), c = (c.coef))
    return(out)
}
print.mKrig <- function(x, digits = 4, ...) {

       if( is.matrix(x$residuals)){
         n<- nrow(x$residuals)
         m<- ncol(x$residuals)}
       else{
         n<- length( x$residuals)
         m<-1 }

    c1 <- "Number of Observations:"
    c2 <- n

    if( m>1){
    c1<- c( c1,"Number of data sets fit:")
    c2<- c( c2, m)}
    c1 <- c(c1, "Degree of polynomial null space ( base model):")
    c2 <- c(c2, x$m - 1)
    c1 <- c(c1, "Number of parameters in the null space")
    c2 <- c(c2, x$nt)
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

    if( m==1){
    c1 <- c(c1, "MLE sigma ")
    c2 <- c(c2, signif(x$shat.MLE, digits))
    c1 <- c(c1, "MLE rho")}

    c2 <- c(c2, signif(x$rho.MLE, digits))
    c1 <- c(c1, "Nonzero entries in covariance")
    c2 <- c(c2, x$nonzero.entries)
    sum <- cbind(c1, c2)
    dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
    cat("Call:\n")
    dput(x$call)
    print(sum, quote = FALSE)
    cat("Covariance Model:", x$cov.function, fill = TRUE)
    if (x$cov.function == "stationary.cov") {
#  this unattractive line is to handle the name of the Covariance
#  in the case it is the default or passed as a function (not quoted).
        cat("  Covariance function is ", ifelse(is.null(x$args$Covariance), 
            "Exponential", x$args$Covariance), fill = TRUE)
    }
    if (!is.null(x$args)) {
        cat("  Names of non-default covariance arguments and their values ", 
            fill = TRUE)
        nlist <- as.character(names(x$args))
        NL <- length(nlist)
        for (k in 1:NL) {
            cat("Argument named: ", nlist[k], " is ", fill=TRUE)
            ifelse(length(x$args[[k]]) < 11, print(x$args[[k]]), 
                cat( "first 10 values: ", 
                      (x$args[[k]])(1:10), "...", fill = TRUE))
        }
    }
    invisible(x)
}

predict.mKrig <- function(object, xnew = NULL, ynew = NULL, 
    derivative = 0, ...) {
    # the main reason to pass new args to the covariance is to increase
    # the temp space size for sparse multiplications
    # other optional arguments from mKrig are passed along in the
    # list object$args
    cov.args <- list(...)
    # predict at observation locations by default
    if (is.null(xnew)) {
        xnew <- object$x
    }
    if (!is.null(ynew)) {
        coef.hold <- mKrig.coef(object, ynew)
        c.coef <- coef.hold$c
        d.coef <- coef.hold$d
    }
    else {
        c.coef <- object$c
        d.coef <- object$d
    }
    # fixed part of the model this a polynomial of degree m-1
    # Tmatrix <- fields.mkpoly(xnew, m=object$m)
    #
    if (derivative == 0) {
        temp1 <- fields.mkpoly(xnew, m = object$m) %*% d.coef
    }
    else {
        temp1 <- fields.derivative.poly(xnew, m = object$m, d.coef)
    }
    # add nonparametric part. Covariance basis functions
    # times coefficients.
    # syntax is the name of the function and then a list with
    # all the arguments. This allows for different covariance functions
    # that have been passed as their name.
    if (derivative == 0) {
        # argument list are the parameters and other options from mKrig
        #  locations and coefficients,
        temp2 <- do.call(object$cov.function.name, c(object$args, 
            list(x1 = xnew, x2 = object$knots, C = c.coef), cov.args))
    }
    else {
        temp2 <- do.call(object$cov.function.name, c(object$args, 
            list(x1 = xnew, x2 = object$knots, C = c.coef, derivative = derivative), 
            cov.args))
    }
    # add two parts together and coerce to vector
    return((temp1 + temp2))
}
