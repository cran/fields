# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
Krig.make.W <- function(out, verbose = FALSE) {
    if (verbose) {
        cat("W", fill = TRUE)
        print(out$W)
    }
    if (out$nondiag.W) {
        #
        # create W from scratch or grab it from passed object
        if (is.null(out$W)) {
            if (verbose) {
                print(out$wght.function.name)
            }
            W <- do.call(out$wght.function.name, c(list(x = out$xM), 
                out$wght.args))
            #       adjust W based on diagonal weight terms
            #
            W <- sqrt(out$weightsM) * t(sqrt(out$weightsM) * 
                W)
        }
        else {
            W <- out$W
        }
        #
        # symmetric square root
        temp <- eigen(W, symmetric = TRUE)
        W2 <- temp$vectors %*% diag(sqrt(temp$values)) %*% t(temp$vectors)
        return(list(W = W, W2 = W2))
    }
    else {
        #
        #  These are created only for use with default method to stay
        #   consistent with nondiagonal elements.
        if (out$fixed.model) {
            return(list(W = NULL, W2 = NULL))
        }
        else {
            return(list(W = out$weightsM, W2 = sqrt(out$weightsM) ))
        }
    }
}
