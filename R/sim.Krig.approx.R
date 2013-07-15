# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"sim.Krig.approx" <- function(object, grid.list = NA, 
    M = 1, nx = 40, ny = 40,  verbose = FALSE, 
     extrap = FALSE) {
    # check that this is a stationary covariance
    if (object$cov.function.name != "stationary.cov") {
        stop("covariance function is not stationary.cov")
    }
    # create grid if not passed
    if (is.na(grid.list)[1]) {
        if (is.null(object$x)) {
            stop("Need a an X matrix in the output object")
        }
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny) 
    }
    #
    # extract what are the x and y and their lengths
    #
    temp <- parse.grid.list(grid.list)
    nx <- temp$nx
    ny <- temp$ny
    #
    # coerce grid list to have x and y components
    #
    glist <- list(x = temp$x, y = temp$y)
    # figure out what sigma and rho should be
   
        sigma2 <- object$best.model[2]
        rho <- object$best.model[3]
    #
    # set up various sizes of arrays
    m <- nx * ny
    n <- nrow(object$xM)
    N <- n
    if (verbose) {
        cat(" m,n,N ", m, n, N, fill = TRUE)
    }
    #transform the new points
    xc <- object$transform$x.center
    xs <- object$transform$x.scale
    #  xpM <- scale(xp, xc, xs)
    if (verbose) {
        cat("center and scale", fill = TRUE)
        print(xc)
        print(xs)
    }
    #
    # set up for simulating on a grid
    #
    cov.obj <- do.call("stationary.image.cov", c(object$args, 
        list(setup = TRUE, grid = glist)))
    out <- array(NA, c(nx, ny, M))
    #
    # find conditional mean field from initial fit
    # don't multiply by sd or add mean if this is
    # a correlation model fit.
    # (these are added at the predict step).
    # from now on all predicted values are on the grid
    # represented by a matrix
    h.hat <- predict.surface(object, grid.list = grid.list, extrap = extrap)$z
    if (verbose) {
        cat("mean predicted field", fill = TRUE)
        image.plot(h.hat)
    }
    # empty surface object to hold ('truth') simulated fields
    h.true <- list(x = glist$x, y = glist$y, z = matrix(NA, nx, 
        ny))
    # covariance matrix for observations
    W2i <- Krig.make.Wi(object, verbose = verbose)$W2i
    if (verbose) {
        cat("dim of W2i", dim(W2i), fill = TRUE)
    }
    ####
    ### begin the big loop
    ###
    for (k in 1:M) {
        # simulate full field
        h.true$z <- sqrt(object$rhohat) * sim.rf(cov.obj)
        if (verbose) {
            cat("mean predicted field", fill = TRUE)
            image.plot(h.true)
        }
        # value of simulated field at observations
        #
        # NOTE: fixed part of model (null space) need not be simulated
        # because the  estimator is unbiased for this part.
        # the variability is still captured because the fixed part
        # is still estimated as part of the predict step below
        #
        #   bilinear interpolation to approximate values at data locations
        #
        h.data <- interp.surface(h.true, object$xM)
        if (verbose) {
            cat("synthetic true values", h.data, fill = TRUE)
        }
        # create synthetic data
        # NOTE:these are actually the 'yM's  the y's
        # having been collapsed to replicate means.
        y.synthetic <- h.data + sqrt(sigma2) * W2i %d*% rnorm(N)
        if (verbose) {
            cat("synthetic data", y.synthetic, fill = TRUE)
        }
        # predict at grid using these data
        # and subtract from 'true' value
        temp.error <- predict.surface(object, grid.list = grid.list, 
            yM = y.synthetic, eval.correlation.model = FALSE, 
            extrap = TRUE)$z - h.true$z
        if (verbose) {
            cat("mean predicted field", fill = TRUE)
            image.plot(temp.error)
        }
        # add the error to the actual estimate  (conditional mean)
        out[, , k] <- h.hat + temp.error
    }
    return(list(x = glist$x, y = glist$y, z = out))
}
