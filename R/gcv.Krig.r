"gcv.Krig" <-
function (out, lambda.grid = NA, cost = 1, nstep.cv = 80, rmse = NA, 
    verbose = FALSE, tol = 1e-05, offset = 0, y = NULL, lambda = NA) 
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
    if (is.na(lambda.grid[1])) {
        l1 <- 1/D[np - nt - 1]
        tr <- np
        for (k in 1:20) {
            tr <- sum(1/(1 + l1 * D))
            if (tr < (nt + 0.05)) 
                break
            l1 <- l1 * 4
        }
        l2 <- 1/max(D)
        for (k in 1:20) {
            tr <- sum(1/(1 + l2 * D))
            cut <- (N - (cost * (tr - nt - offset) + nt) * (N/length(D)))
            if ((tr > (np - offset) * 0.95) | (cut/N <= 0.05)) 
                break
            l2 <- l2/4
        }
        lambda.grid <- exp(seq(log(l2), log(l1), , nstep.cv))
    }
    nl <- length(lambda.grid)
    nd <- length(D)
    trA <- MSE <- RSS.model <- rep(NA, nl)
    Dl <- rep(NA, nd)
    for (k in 1:nl) {
        Dl <- lambda.grid[k] * D
        RSS.model[k] <- sum(((u * Dl)/(1 + Dl))^2)
        trA[k] <- sum(1/(1 + Dl))
    }
    RSS <- pure.ss + RSS.model
    MSE <- pure.ss/N + RSS.model/length(D)
    MSE.one <- (pure.ss + RSS.model)/N
    den <- 1 - (cost * (trA - nt - offset) + nt)/length(D)
    den.one <- 1 - (cost * (trA - nt - offset) + nt)/N
    V <- ifelse(den > 0, (MSE)/den^2, NA)
    V.one <- ifelse(den.one > 0, MSE.one/den.one^2, NA)
    V.model <- ifelse(den > 0, ((RSS.model/length(D))/den^2), 
        NA)
    shat <- sqrt(RSS/(N - trA))
    gcv.grid <- cbind(lambda.grid, trA, V, V.one, V.model, shat)
    gcv.grid <- as.data.frame(gcv.grid)
    names(gcv.grid) <- c("lambda", "trA", "GCV", "GCV.one", "GCV.model", 
        "shat")
    info <- list(matrices = list(D = D, u = u), N = N, nt = nt, 
        cost = cost, pure.ss = pure.ss, shat.pure.error = shat.pure.error, 
        offset = offset)
    if (verbose) 
        print(info)
    lambda.est <- matrix(ncol = 4, nrow = 5, dimnames = list(c("GCV", 
        "GCV.model", "GCV.one", "RMSE", "pure error"), c("lambda", 
        "trA", "GCV", "shat")))
    lambda.est[1, 1] <- Krig.find.gcvmin(info, lambda.grid, gcv.grid$GCV, 
        Krig.fgcv, tol = tol, verbose = verbose)
    if (!is.na(shat.pure.error)) {
        temp <- gcv.grid$GCV.model
        lambda.est[2, 1] <- Krig.find.gcvmin(info, lambda.grid, 
            temp, Krig.fgcv.model, tol = tol, verbose = verbose)
    }
    lambda.est[3, 1] <- Krig.find.gcvmin(info, lambda.grid, gcv.grid$GCV.one, 
        Krig.fgcv.one, tol = tol, verbose = verbose)
    lambda.rmse <- NA
    lambda.pure.error <- NA
    if (!is.na(rmse)) {
        if( all(gcv.grid$shat < rmse)| all(gcv.grid$shat >rmse)) {guess<- NA}
        else{        guess <- max(gcv.grid$lambda[gcv.grid$shat < rmse])}
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
            warning("Value of rmse is outside possible range")
        }
    }
    if (!is.na(shat.pure.error)) {
#        guess <- max(gcv.grid$lambda[gcv.grid$shat < shat.pure.error])
if( all(gcv.grid$shat < shat.pure.error)| all(gcv.grid$shat >shat.pure.error)) 
             {guess<- NA}
        else{ 
             guess <- max(gcv.grid$lambda[gcv.grid$shat < shat.pure.error])}

        if (!is.na(guess)&(guess!= -Inf)) {
            lambda.pure.error <- find.upcross(Krig.fs2hat, info, 
                upcross.level = shat.pure.error^2, guess = guess, 
                tol = tol * shat.pure.error^2)
            lambda.est[5, 1] <- lambda.pure.error
        }
        else {
            warning("Value of pure error estimate  is outside possible range")
        }
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
    list(gcv.grid = gcv.grid, lambda.est = lambda.est, lambda.best = lambda.est[1, 
        1])
}
