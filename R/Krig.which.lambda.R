Krig.which.lambda <- function(out) {
    #
    # determine the method for finding lambda
    #  Note order
    # default is to do 'gcv/REML'
    out2 <- list()
    # copy all all parameters to out2 just to make this
    # easier to read.
    out2$method <- out$method
    out2$lambda.est <- NA
    out2$lambda <- out$lambda
    out2$eff.df <- out$eff.df
    out2$rho <- out$rho
    out2$sigma2 <- out$sigma2
    if (!is.na(out2$lambda) | !is.na(out2$eff.df)) {
        #
        # this indicates lambda has been supplied and leads to
        # the cholesky type computational approaches
        #        -- but only if GCV is FALSE
        #
        out2$method <- "user"
    }
    out2$GCV <- out$GCV
    if (!is.na(out2$eff.df)) {
        #
        # this indicates df has been supplied and needs
        # GCV to be true to compute the lambda
        # that matches the df
        #
        out2$GCV <- TRUE
    }
    if (!is.na(out2$rho) & !is.na(out2$sigma2)) {
        out2$method <- "user"
        out2$lambda <- out2$sigma2/out2$rho
    }
    #
    # NOTE: method='user' means that a value of lambda has been supplied
    #        and so GCV etc to determine lambda is not needed.
    #  gcv TRUE means that the decompositions will be done to
    #    evaluate the estimate at arbitrary lambda (and also be
    #    able to compute the effective degrees of freedom).
    #
    #    The fixed lambda calculations are very efficient but
    #    do not make it feasible for GCV/REML  or effective degrees of
    #    freedom calculations.
    #
    out2$fixed.model <- (out2$method == "user") & (!out2$GCV)
    #
    return(out2)
}
