# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"Krig.make.u" <- function(out, y = NULL, yM = NULL, 
    verbose = FALSE) {
    #
    # Determine whether to collapse onto means of replicates ( using y)
    # if the data has been passed use as the replicate means (yM) use that.
    # If both y and YM are null then just use out$yM
    # For readability of this function, all this tortured logic happens in
    #  Krig.ynew.
    #
    out2 <- Krig.ynew(out, y, yM)
    temp.yM <- out2$yM
    nt <- out$nt
    np <- out$np
    ndata <- ncol(temp.yM)
    u <- NA
    call.name <- out$cov.function.name
    if (verbose) {
        cat("dimension of yM in Krig.coef", fill = TRUE)
        print(dim(temp.yM))
    }
    #
    #   case when knots= unqiue x's
    # any lambda
    #
    if (out$decomp == "WBW") {
        # pad u with zeroes that corresond to null space basis functions
        # this makes it compatible with the DR decomposition.
        u <- rbind(matrix(0, nrow = out$nt, ncol = ndata), t(out$matrices$V) %*% 
            qr.q2ty(out$matrices$qr.T, out$W2 %d*% temp.yM))
    }
    #
    # case with knots
    # any lambda
    #
    if (out$decomp == "DR") {
        # X is the monster matrix ...  X = [ M | K]
        X <- cbind(do.call(out$null.function.name, c(out$null.args, 
            list(x = out$xM, Z = out$ZM))), do.call(call.name, 
            c(out$args, list(x1 = out$xM, x2 = out$knots))))
        u <- t(out$matrices$G) %*% t(X) %*% (out$weightsM %d*% 
            temp.yM)
    }
    return(list(u = u, shat.rep = out2$shat.rep, shat.pure.error = out2$shat.pure.error, 
        pure.ss = out2$pure.ss))
}
