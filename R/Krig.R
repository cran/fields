# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"Krig" <- function(x, Y, cov.function = "stationary.cov", 
    lambda = NA, df = NA, GCV = FALSE, Z = NULL, cost = 1, knots = NA, 
    weights = NULL, m = 2, nstep.cv = 200, scale.type = "user", 
    x.center = rep(0, ncol(x)), x.scale = rep(1, ncol(x)), rho = NA, 
    sigma2 = NA, method = "GCV", verbose = FALSE, mean.obj = NA, 
    sd.obj = NA, null.function = "Krig.null.function", wght.function = NULL, 
    offset = 0, outputcall = NULL, na.rm = TRUE, cov.args = NULL, 
    chol.args = NULL, null.args = NULL, wght.args = NULL, W = NULL, 
    give.warnings = TRUE, ...) # NOTES
# the verbose switch prints many intermediate steps as an aid in debugging.
#
{
    #
    # create output list
    out <- list()
    ###########################################################
    #  First series of steps simply store pieces of the passed
    #    information to output list (i.e. the Krig object)
    ##########################################################
    if (is.null(outputcall)) {
        out$call <- match.call()
    }
    else {
        out$call <- outputcall
    }
    #
    # save covariance function as its name
    #
    out$cov.function.name <- as.character(substitute(cov.function))
    #
    # save null space function as its name
    #
    out$null.function.name <- as.character(substitute(null.function))
    #
    # save weight  function as its name if it is not a NULL
    #
    if (is.null(wght.function)) {
        out$wght.function.name <- NULL
    }
    else {
        out$wght.function.name <- as.character(substitute(wght.function))
    }
    out$W <- W
    if (verbose) {
        print(out$cov.function.name)
        print(out$null.function.name)
        print(out$wght.function.name)
    }
    #
    # logical to indicate if the 'C' argument is present in cov.function
    #
    C.arg.missing <- all(names(formals(get(out$cov.function.name))) != 
        "C")
    if (C.arg.missing) 
        stop("Need to have C argument in covariance function\nsee Exp.cov.simple as an example")
    #
    # save parameters values possibly passed to the covariance function
    # also those added to call are assumed to be covariance arguments.
    if (!is.null(cov.args)) 
        out$args <- c(cov.args, list(...))
    else out$args <- list(...)
    #
    # default values for null space function
    out$null.args <- null.args
    #
    #       set degree of polynomial null space if this is default
    #       mkpoly is used so often is it helpful to include m argument
    #       by default in Krig call.
    if (out$null.function.name == "Krig.null.function") {
        out$null.args <- list(m = m)
        out$m <- m
    }
    #
    # default values for Cholesky decomposition, these are important
    # for sparse matrix decompositions used in Krig.engine.fixed.
    if (is.null(chol.args)) {
        out$chol.args <- list(pivot = FALSE)
    }
    else {
        out$chol.args <- chol.args
    }
    # additional arguments for weight matrix.
    out$wght.args <- wght.args
    #
    # the offset is the effective number of parameters used in the GCV
    # calculations -- unless this is part of an additive model this
    # is likely zero
    out$offset <- offset
    #
    # the cost is the multiplier applied to the GCV eff.df
    # lambda and df are two ways of parameterizing the smoothness
    # and are related by a monotonic function that unfortunately
    # depends on the locations of the data.
    # lambda can be used directly in the linear algebra, df
    # must be transformed to lambda numerically using the monotonic trransformation
    # sigma2 is the error variance and rho the multiplier for the covariance
    # method is how to determine lambda
    # the GCV logical forces the code to do the more elaborate decompositions
    # that faclitate estimating lambda -- even if a specific lambda value is
    # given.
    out$cost <- cost
    out$lambda <- lambda
    out$eff.df <- df
    out$sigma2 <- sigma2
    out$rho <- rho
    out$method <- method
    out$GCV <- GCV
    #
    # correlation model information
    #
    out$mean.obj <- mean.obj
    out$sd.obj <- sd.obj
    out$correlation.model <- !(is.na(mean.obj[1]) & is.na(sd.obj[1]))
    #
    # transformation info
    out$scale.type <- scale.type
    out$x.center <- x.center
    out$x.scale <- x.scale
    #
    # verbose block
    if (verbose) {
        cat("  Cov function arguments in call  ", fill = TRUE)
        print(out$args)
        cat(" covariance function used is : ", fill = TRUE)
        print(out$cov.function.name)
    }
    ###############################################################
    # Begin modifications and transformations of input information
    # note that many of these manipulations follow a strategy
    # of passing the Krig object (out) to a function and
    # then appending the information from this function to
    # the Krig object. In this way the Krig object  is built up
    # in steps and it is hoped easier to follow.
    ###############################################################
    # various checks on x and  Y including removal of NAs in Y
    if (verbose) {
        cat("checks on x,Y, and Z", fill = TRUE)
    }
    # Here is an instance of adding to the Krig object
    # in this case some onerous bookkeeping making sure arguments are consistent
    out2 <- Krig.check.xY(x, Y, Z, weights, na.rm, verbose = verbose)
    out <- c(out, out2)
    # transform to correlation model (if appropriate)
    # find replicates and collapse to means and pool variances.
    # Transform unique x locations and knots.
    if (out$correlation.model) {
        out$y <- Krig.cor.Y(out, verbose = verbose)
    }
    if (verbose) {
        cat("Transform x and knots: ", fill = TRUE)
    }
    out2 <- Krig.transform.xY(out, knots, verbose = verbose)
    out <- c(out, out2)
    # NOTE: knots have been transformed after this step
    #############################################################
    #  Figure out what to do
    #############################################################
    #
    # this functions works through the logic of
    # what has been supplied for lambda
    out2 <- Krig.which.lambda(out)
    out[names(out2)] <- out2
    # verbose block
    if (verbose) {
        cat("Modified values for smoothing controls", fill = TRUE)
        print(out2)
    }
    # verbose block
    if (verbose) {
        cat("lambda, fixed? ", lambda, out$fixed.model, fill = TRUE)
    }
    # Make weight matrix for observations
    #    ( this is proportional to the inverse square root of obs covariance)
    #     if a weight function or W has not been passed then this is
    #     diag( out$weightsM) for W
    #     The checks represent a limitation of this model to
    #     the  WBW type decoposition and no replicate observations.
    out$nondiag.W <- (!is.null(wght.function)) | (!is.null(W))
    if (verbose) {
        cat("out$nondiag", out$nondiag, fill = TRUE)
    }
    # Do not continue if there there is a nondiagonal weight matrix
    # and replicate observations.
    if (out$nondiag.W) {
        if (out$knot.model | out$fixed.model) {
            stop("Non diagonal weight matrix for observations not supported\nwith knots or fixed lambda.")
        }
        if (!is.na(out$shat.pure.error)) {
            stop("Non diagonal weight matrix not implemented with replicate\nlocations")
        }
    }
    #  make weight matrix and its square root having passed checks
    out <- c(out, Krig.make.W(out, verbose = verbose))
    ########################################################
    #  You have reached the Engines!
    ########################################################
    #   Do the intensive linear algebra to find the solutions
    #   this is where all the heavy lifting happens.
    #
    #   Note that all the information is passed as a list
    #   including arguments to the cholesky decomposition
    #   used within Krig.engine.fixed
    #
    # The results are saved in the component matrices
    #
    # if method=='user' then just evaluate at single lambda
    #  fixed here means a fixed lambda
    #
    # For fixed lambda the decompositions with and without knots
    # are surprisingly similar and so are in one engine.
    ###########################################################
    if (out$fixed.model) {
        out$matrices <- Krig.engine.fixed(out, verbose = verbose)
        # can't find the trace of A matrix in fixed lambda case so set this to NA.
        out$eff.df <- NA
    }
    #
    # alternative are
    # matrix decompositions suitable for
    # evaluation at many lambdas to facilitate GCV/REML estimates  etc.
    #
    if (!out$fixed.model) {
        if (out$knot.model) {
            # the knot model engine
            out$matrices <- Krig.engine.knots(out, verbose = verbose)
            out$pure.ss <- out$matrices$pure.ss
        }
        else {
            # standard engine following the basic computations for thin plate splines
            if (verbose) {
                cat("Call to Krig.engine.default", fill = TRUE)
            }
            out$matrices <- Krig.engine.default(out, verbose = verbose)
        }
    }
    #
    # store basic information about decompositions
    out$nt <- out$matrices$nt
    out$np <- out$matrices$np
    out$decomp <- out$matrices$decomp
    #
    # Now determine a logical vector indices for coefficients tied to  the
    # the 'spatial drift' i.e. the fixed part of the model
    # that is not due to the Z covariates.
    # NOTE that the spatial drift coefficients must be the first columns of the
    # M matrix
    if (is.null(out$Z)) {
        out$ind.drift <- rep(TRUE, out$nt)
    }
    else {
        mZ <- ncol(out$ZM)
        out$ind.drift <- c(rep(TRUE, out$nt - mZ), rep(FALSE, 
            mZ))
    }
    if (verbose) {
        cat("null df: ", out$nt, "drift df: ", sum(out$ind.drift), 
            fill = TRUE)
    }
    #########################
    # End of engine block
    #########################
    #################################################
    # Do GCV and REML search over lambda if not fixed
    #################################################
    if (!out$fixed.model | out$GCV) {
        if (verbose) {
            cat("call to gcv.Krig", fill = TRUE)
        }
        gcv.out <- gcv.Krig(out, nstep.cv = nstep.cv, verbose = verbose, 
            cost = out$cost, offset = out$offset, give.warnings = give.warnings)
        out$gcv.grid <- gcv.out$gcv.grid
        #
        #  a handy summary table of the search results
        out$lambda.est <- gcv.out$lambda.est
        #
        # verbose block
        if (verbose) {
            cat("returned GCV and REML grid search", fill = TRUE)
            print(out$gcv.grid)
        }
        #
        # assign the preferred lambda either from GCV/REML/MSE or the user value
        # NOTE: gcv/reml can be done but the estimate is
        # still evaluted at the passed  user values of lambda (or df)
        # If df is passed need to calculate the implied lambda value
        if (out$method != "user") {
            out$lambda <- gcv.out$lambda.est[out$method, 1]
            out$eff.df <- out$lambda.est[out$method, 2]
        }
        else {
            if (!is.na(out$eff.df)) {
                out$lambda <- Krig.df.to.lambda(out$eff.df, out$matrices$D)
            }
            else {
                out$eff.df <- Krig.ftrace(out$lambda, out$matrices$D)
            }
        }
        if (verbose) {
            cat("lambda set in GCV block", fill = TRUE)
            print(out$lambda)
            cat("trace of A", fill = TRUE)
            print(out$eff.df)
        }
    }
    ##########################
    # end GCV/REML block
    ##########################
    #
    # Now we clean up what has happen and stuff into output object.
    #
    ##########################################
    # find coefficients at prefered lambda
    # and evaluate the solution at observations
    ##########################################
    #   pass replicate group means -- no need to recalculate these.
    if (verbose) {
        cat("Call to Krig.coef:", fill = TRUE)
    }
    out2 <- Krig.coef(out, yM = out$yM, verbose = verbose)
    out <- c(out, out2)
    if (verbose) {
        cat("Krig.coef:", fill = TRUE)
        print(out2)
    }
    #######################################################################
    # fitted values and residuals and predicted values for full model and
    # also on the null space (fixed
    # effects). But be sure to do this at the nonmissing x's.
    ##################################################################
    out$fitted.values <- predict.Krig(out, x = out$x, Z = out$Z, 
        eval.correlation.model = FALSE)
    out$residuals <- out$y - out$fitted.values
    #
    # this is just M%*%d  note use of do.call using function name
    Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
        list(x = out$x, Z = out$Z)))
    out$fitted.values.null <- as.matrix(Tmatrix) %*% out$d
    #
    # verbose block
    if (verbose) {
        cat("residuals", out$residuals, fill = TRUE)
    }
    #
    # find various estimates of sigma and rho
    out2 <- Krig.parameters(out)
    out <- c(out, out2)
    ################################################
    # assign the 'best' model as a default choice
    # either use the user supplied values or the results from
    # optimization
    ################################################
    passed.sigma2 <- (!is.na(out$sigma2))
    if (out$method == "user" & passed.sigma2) {
        out$best.model <- c(out$lambda, out$sigma2, out$rho)
    }
    else {
        # in this case lambda is from opt. or supplied by user
        out$best.model <- c(out$lambda, out$shat.MLE^2, out$rhohat)
    }
    # Note: values in best.model are used in subsquent functions as the choice
    # for these parameters!
    # set class
    class(out) <- c("Krig")
    return(out)
}
