"gcv.Krig" <-
function (out, lambda.grid = NA, cost = 1, nstep.cv = 80, rmse = NA, 
    verbose = FALSE, tol = 1e-05, offset = 0, y = NULL, 
give.warnings = TRUE,give.warnings.REML = FALSE) 
{
    nt <- out$nt
    np <- out$np
    N <- out$N
    D <- out$matrices$D
    if (is.null(y)) {
        u <- out$matrices$u
        shat.pure.error <- out$shat.pure.error
        pure.ss <- out$pure.ss
    }
    else {
        out2 <- Krig.updateY(out, y)
        u <- out2$u
        shat.pure.error <- out2$shat.pure.error
        pure.ss <- out2$pure.ss
    }
# 
# set limits on serach for lambda 
# this is best done realtive to the trace (effective degrees of freeedom) 
#
   if (is.na(lambda.grid[1])) {
        temp.df <- seq(nt, (np - offset) * 0.95, , nstep.cv)
        temp.df[1] <- temp.df[1] + 0.001
        for (k in 1:nstep.cv) {
            lambda.grid[k] <- Krig.df.to.lambda(temp.df[k], D)
        }
    }
    lambda.grid <- sort(lambda.grid)
    nl <- length(lambda.grid)
    nd <- length(D)
    V <- V.model <- V.one <- lplike <- trA <- shat <- rep(NA, 
        nl)
    Dl <- rep(NA, nd)
    info <- list(matrices = list(D = D, u = u), N = N, nt = nt, 
        cost = cost, pure.ss = pure.ss, shat.pure.error = shat.pure.error, 
        offset = offset)
#
# grid search to get rough idea of lambda estimates. 
# 
   for (k in 1:nl) {
        V[k] <- Krig.fgcv(lambda.grid[k], info)
        V.one[k] <- Krig.fgcv.one(lambda.grid[k], info)
        V.model[k] <- Krig.fgcv.model(lambda.grid[k], info)
        lplike[k] <- Krig.flplike(lambda.grid[k], info)
        shat[k] <- sqrt(Krig.fs2hat(lambda.grid[k], info))
        trA[k] <- Krig.ftrace(lambda.grid[k], D)
    }
#
# save results in matrix
#
    gcv.grid <- cbind(lambda.grid, trA, V, V.one, V.model, shat, 
        lplike)
    gcv.grid <- as.data.frame(gcv.grid)
    names(gcv.grid) <- c("lambda", "trA", "GCV", "GCV.one", "GCV.model", 
        "shat", "-Log Profile ")
    if (verbose) 
        print(info)
    if (verbose) 
        print(gcv.grid)
    lambda.est <- matrix(ncol = 4, nrow = 6, dimnames = list(c("GCV", 
        "GCV.model", "GCV.one", "RMSE", "pure error", "-Log Profile (REML)"), 
        c("lambda", "trA", "GCV", "shat")))
#
# find various estimates for the smoothing parameter
# using grid serach as starting values. 
#
    lambda.est[1, 1] <- Krig.find.gcvmin(info, lambda.grid, gcv.grid$GCV, 
                        Krig.fgcv, tol = tol, verbose = verbose,
                        give.warnings= give.warnings)
    if (!is.na(shat.pure.error)) {
        temp <- gcv.grid$GCV.model
        lambda.est[2, 1] <- Krig.find.gcvmin(info, lambda.grid, 
                            temp, Krig.fgcv.model, tol = tol, verbose = verbose,
                            give.warnings= give.warnings)
    }
    lambda.est[3, 1] <- Krig.find.gcvmin(info, lambda.grid, gcv.grid$GCV.one, 
                        Krig.fgcv.one, tol = tol, verbose = verbose,
                        give.warnings= give.warnings)

    lambda.est[6, 1] <- Krig.find.REML(info, lambda.grid, 
                        gcv.grid$"-Log Profile", 
                        Krig.flplike, tol = tol, verbose = verbose,
                        give.warnings= give.warnings.REML)
    if (verbose) {
        cat(" mle estimate", lambda.est[6, 1], fill = TRUE)
    }
    lambda.rmse <- NA
    lambda.pure.error <- NA
    if (!is.na(rmse)) {
        if (all(gcv.grid$shat < rmse) | all(gcv.grid$shat > rmse)) {
            guess <- NA
        }
        else {
            guess <- max(gcv.grid$lambda[gcv.grid$shat < rmse])
        }
        if (verbose) {
            print(rmse)
            print(guess)
        }
        if (!is.na(guess)) {
            lambda.rmse <- find.upcross(Krig.fs2hat, info, upcross.level = rmse^2, 
                guess = guess, tol = tol * rmse^2)
            lambda.est[4, 1] <- lambda.rmse
        }
        else {
            if (give.warnings) {
                warning("Value of rmse is outside possible range")
            }
        }
    }
    if (!is.na(shat.pure.error)) {
        if (all(gcv.grid$shat < shat.pure.error) | all(gcv.grid$shat > 
            shat.pure.error)) {
            guess <- NA
        }
        else {
            guess <- max(gcv.grid$lambda[gcv.grid$shat < shat.pure.error])
        }
        if (!is.na(guess) & (guess != -Inf)) {
            lambda.pure.error <- find.upcross(Krig.fs2hat, info, 
                upcross.level = shat.pure.error^2, guess = guess, 
                tol = tol * shat.pure.error^2)
            lambda.est[5, 1] <- lambda.pure.error
        }
        else {
            if (give.warnings) {
                warning("Value of pure error estimate  is outside possible range")
            }
        }
    }
    if (verbose) {
        cat(" Done with refining", fill = TRUE)
    }
    for (k in 1:5) {
        lam <- lambda.est[k, 1]
        if (!is.na(lam)) {
            lambda.est[k, 2] <- Krig.ftrace(lam, D)
            if (k == 1 | k > 3) {
                lambda.est[k, 3] <- Krig.fgcv(lam, info)
            }
            if (k == 2) {
                lambda.est[k, 3] <- Krig.fgcv.model(lam, info)
            }
            if (k == 3) {
                lambda.est[k, 3] <- Krig.fgcv.one(lam, info)
            }
            lambda.est[k, 4] <- sqrt(Krig.fs2hat(lam, info))
        }
    }
    lam.ml <- lambda.est[6, 1]
    lambda.est[6, 2] <- Krig.ftrace(lam.ml, D)
    lambda.est[6, 3] <- Krig.fgcv(lam.ml, info)
    lambda.est[6, 4] <- sqrt(Krig.fs2hat(lam.ml, info))
    list(gcv.grid = gcv.grid, lambda.est = lambda.est, lambda.best = lambda.est[1, 
        1])
}

