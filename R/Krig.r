"Krig" <-
function (x, Y, cov.function = exp.cov, lambda = NA, df = NA, 
    cost = 1, knots, weights = rep(1, length(Y)), m = 2, return.matrices = T, 
    nstep.cv = 80, scale.type = "user", x.center = rep(0, ncol(x)), 
    x.scale = rep(1, ncol(x)), rho = NA, sigma2 = NA, method = "GCV", 
    decomp = "DR", verbose = F, cond.number = 10^8, mean.obj = NULL, 
    sd.obj = NULL, yname = NULL, return.X = T, null.function = make.tmatrix, 
    offset = 0, outputcall = NULL, cov.by.name = T, ...) 
{
    out <- list()
    class(out) <- c("Krig")
    if (is.null(outputcall)) {
        out$call <- match.call()
    }
    else {
        out$call <- outputcall
    }
    N <- length(Y)
    if (sum(is.na(Y)) != 0) {
        stop("Need to remove missing values\n from Y vector!")
    }
    out$cov.by.name <- cov.by.name
    if (is.null(yname)) 
        out$yname <- deparse(substitute(Y))
    else out$yname <- yname
    out$N <- N
    out$decomp <- decomp
    out$make.tmatrix <- null.function
    out$offset <- offset
    out$cost <- cost
    out$call.name <- as.character(substitute(cov.function))
    if (!cov.by.name) {
        out$cov.function <- replace.args.function(cov.function, 
            ...)
    }
    if (cov.by.name) {
        out$args <- list(...)
        if (verbose) {
            cat(" Local cov function arguments ", fill = T)
            print(out$args)
            cat(" covariance function used is : ", fill = T)
            print(out$call.name)
        }
    }
    if (verbose & !cov.by.name) {
        cat(" Local cov function arguments", fill = T)
        print(args(out$cov.function))
    }
    x <- as.matrix(x)
    if (length(dimnames(x)) != 2) {
        dimnames(x) <- list(format(1:nrow(x)), paste("X", format(1:ncol(x)), 
            sep = ""))
    }
    if (length(dimnames(x)[[1]]) == 0) {
        dimnames(x)[[1]] <- format(1:nrow(x))
    }
    if (length(dimnames(x)[[2]]) == 0) {
        dimnames(x)[[2]] <- paste("X", format(1:ncol(x)), sep = "")
    }
    Y <- c(Y)
    out$y <- Y
    out$x <- x
    if (!is.null(sd.obj) & !is.null(mean.obj)) {
        correlation.model <- T
    }
    else correlation.model <- F
    out$correlation.model <- correlation.model
    if (correlation.model) {
        Yraw <- Y
        out$mean.obj <- mean.obj
        out$sd.obj <- sd.obj
        Y <- (Y - predict(mean.obj, x))/predict(sd.obj, x)
        if (verbose) 
            print(Y)
    }
    out$weights <- weights
    rep.info <- cat.matrix(x)
    out$rep.info <- rep.info
    if (verbose) {
        cat("rep.info", fill = T)
        print(rep.info)
    }
    if (max(rep.info) == N) {
        shat.rep <- NA
        shat.pure.error <- NA
        out$pure.ss <- 0
        YM <- Y
        weightsM <- weights
        xM <- as.matrix(x[!dup(rep.info), ])
    }
    else {
        rep.info.aov <- fast.1way(rep.info, Y, weights)
        shat.pure.error <- sqrt(rep.info.aov$MSE)
        shat.rep <- shat.pure.error
        YM <- rep.info.aov$means
        weightsM <- rep.info.aov$w.means
        xM <- as.matrix(x[!dup(rep.info), ])
        out$pure.ss <- rep.info.aov$SSE
        if (verbose) {
            cat(" rep info", fill = T)
            print(rep.info.aov)
        }
    }
    out$yM <- YM
    out$xM <- xM
    out$weightsM <- weightsM
    out$shat.rep <- shat.rep
    out$shat.pure.error <- shat.pure.error
    if (missing(knots)) {
        knots <- xM
        mle.calc <- T
        knot.model <- F
    }
    else {
        mle.calc <- F
        knot.model <- T
    }
    out$knot.model <- knot.model
    out$mle.calc <- mle.calc
    knots <- as.matrix(knots)
    out$knots <- knots
    xM <- transformx(xM, scale.type, x.center, x.scale)
    transform <- attributes(xM)
    if (verbose) {
        cat("transform", fill = T)
        print(transform)
    }
    knots <- scale(knots, center = transform$x.center, scale = transform$x.scale)
    if (verbose) {
        cat("knots in transformed scale", fill = T)
        print(knots)
    }
    out$transform <- transform
    if (!is.na(lambda) | !is.na(df)) 
        method <- "user"
    if (!is.na(rho) & !is.na(sigma2)) {
        lambda <- sigma2/rho
        method <- "user"
    }
    just.solve <- (lambda[1] == 0)
    if (is.na(just.solve)) 
        just.solve <- F
    if (verbose) 
        cat("lambda", lambda, fill = T)
    d <- ncol(xM)
    if (decomp == "DR") {
        qr.T <- qr(out$make.tmatrix(knots, m))
        if (verbose) {
            print(qr.T)
        }
        if (!out$cov.by.name) {
            tempM <- qr.yq2(qr.T, out$cov.function(knots, knots))
        }
        if (out$cov.by.name) {
            tempM <- qr.yq2(qr.T, do.call(out$call.name, c(out$args, 
                list(x1 = knots, x2 = knots))))
        }
        tempM <- qr.q2ty(qr.T, tempM)
        if (verbose) {
            print(dim(tempM))
        }
    }
    if (decomp == "WBW") {
        qr.T <- qr(sqrt(weightsM) * out$make.tmatrix(knots, m))
        if (verbose) {
            print(qr.T)
        }
        if (!out$cov.by.name) {
            tempM <- sqrt(weightsM) * t(sqrt(weightsM) * t(out$cov.function(knots, 
                knots)))
        }
        if (out$cov.by.name) {
            tempM <- sqrt(weightsM) * t(sqrt(weightsM) * t(do.call(out$call.name, 
                c(out$args, list(x1 = knots, x2 = knots)))))
        }
        tempM <- qr.yq2(qr.T, tempM)
        tempM <- qr.q2ty(qr.T, tempM)
    }
    np <- nrow(knots)
    nt <- (qr.T$rank)
    out$np <- np
    out$nt <- nt
    if (verbose) 
        cat("np, nt", np, nt, fill = T)
    if (verbose) 
        print(knots)
    if (just.solve) {
        if (!out$cov.by.name) {
            beta <- qr.coef(qr(cbind(out$make.tmatrix(xM, m), 
                qr.yq2(qr.T, out$cov.function(xM, knots)))), 
                YM)
        }
        if (out$cov.by.name) {
            beta <- qr.coef(qr(cbind(out$make.tmatrix(xM, m), 
                qr.yq2(qr.T, do.call(out$call.name, c(out$args, 
                  list(x1 = xM, x2 = knots)))))), YM)
        }
    }
    else {
        if (decomp == "DR") {
            if (verbose) 
                cat("Type of decomposition", decomp, fill = T)
            H <- matrix(0, ncol = np, nrow = np)
            H[(nt + 1):np, (nt + 1):np] <- tempM
            if (!out$cov.by.name) {
                X <- cbind(out$make.tmatrix(xM, m), qr.yq2(qr.T, 
                  out$cov.function(xM, knots)))
            }
            if (out$cov.by.name) {
                X <- cbind(out$make.tmatrix(xM, m), qr.yq2(qr.T, 
                  do.call(out$call.name, c(out$args, list(x1 = xM, 
                    x2 = knots)))))
            }
            if (verbose) {
                print(weightsM)
                cat("first 3 rows of X", fill = T)
                print(X[1:3, ])
            }
            temp <- svd(sqrt(weightsM) * X)[c("v", "d")]
            cond.matrix <- max(temp$d)/min(temp$d)
            if (cond.matrix > cond.number) {
                print(paste("condition number is ", cond.matrix, 
                  sep = ""))
                stop("Covariance matrix is close\nto\nsingular")
            }
            B <- temp$v %*% diag(1/(temp$d)) %*% t(temp$v)
            temp <- svd(B %*% H %*% B)
            U <- temp$u
            D <- temp$d
            if (verbose) {
                cat("singular values:", fill = T)
                print(D)
            }
            D[(1:nt) + (np - nt)] <- 0
            G <- B %*% U
            u <- t(G) %*% t(X) %*% (weightsM * YM)
            if (verbose) {
                cat("DR: u", fill = T)
                print(u)
            }
            out$pure.ss <- sum(weightsM * (YM - X %*% G %*% u)^2) + 
                out$pure.ss
            if (verbose) {
                cat("pure.ss", fill = T)
                print(out$pure.ss)
            }
            out$matrices <- list(B = B, U = U, u = u, D = D, 
                G = G, qr.T = qr.T, X = X)
        }
        if (decomp == "WBW") {
            temp <- svd(tempM)[c("d", "v")]
            D <- c(rep(0, nt), 1/temp$d)
            if (verbose) {
                cat("singular values:", fill = T)
                print(D)
            }
            G <- matrix(0, ncol = np, nrow = np)
            G[(nt + 1):np, (nt + 1):np] <- temp$v
            G <- G * matrix(D, ncol = np, nrow = np, byrow = T)
            u <- c(rep(0, nt), t(temp$v) %*% qr.q2ty(qr.T, sqrt(weightsM) * 
                YM))
            if (verbose) {
                cat("WBW: u", fill = T)
                print(u)
            }
            if (verbose) {
                cat("WBW: pure.ss", fill = T)
                print(out$pure.ss)
            }
            out$matrices <- list(u = u, D = D, G = G, qr.T = qr.T, 
                decomp = decomp, V = temp$v)
        }
        if (verbose) {
            cat("call to gcv minimization", fill = T)
        }
        gcv.out <- gcv.Krig(out, nstep.cv = nstep.cv, verbose = verbose, 
            cost = out$cost, offset = out$offset)
        gcv.grid <- gcv.out$gcv.grid
        out$gcv.grid <- gcv.grid
        out$lambda.est <- gcv.out$lambda.est
        if (verbose) {
            print(out$gcv.grid)
        }
    }
    if (method == "user") {
        if (!is.na(df)) 
            lambda <- Krig.df.to.lambda(df, D)
        lambda.best <- lambda
    }
    else {
        lambda.best <- gcv.out$lambda.est[method, 1]
    }
    beta <- G %*% ((1/(1 + lambda.best * D)) * u)
    out$m <- m
    if (!just.solve) {
        out$eff.df <- sum(1/(1 + lambda.best * D))
        out$trace <- out$eff.df
        if (verbose) {
            cat("trace of A", fill = T)
            print(out$trace)
        }
    }
    else {
        out$eff.df <- out$np
    }
    if (just.solve) 
        out$lambda <- lambda
    else out$lambda <- lambda.best
    out$beta <- beta
    hold <- Krig.coef(out)
    out$c <- hold$c
    out$d <- hold$d
    if (verbose) {
        cat(names(out), fill = T)
    }
    out$fitted.values <- predict.Krig(out, x, eval.correlation.model = F)
    out$residuals <- Y - out$fitted.values
    if (verbose) {
        cat("resdiuals", out$residuals, fill = T)
    }
    out$fitted.values.null <- as.matrix(out$make.tmatrix(x, m)) %*% 
        out$d
    out$just.solve <- just.solve
    out$shat.GCV <- sqrt(sum(out$weights * out$residuals^2)/(length(Y) - 
        out$trace))
    if (mle.calc) {
        out$rhohat <- sum(out$c * out$yM)/(N - nt)
        if (is.na(rho)) {
            out$rho <- out$rhohat
        }
        else {
            out$rho <- rho
        }
        if (is.na(sigma2)) 
            sigma2 <- out$rho * out$lambda
        out$sigma2 <- sigma2
        out$shat.MLE <- sqrt(out$rhohat * out$lambda)
        out$best.model <- c(out$lambda, out$sigma2, out$rho)
    }
    if (!mle.calc) {
        out$sigma2 <- out$shat.GCV^2
        out$rho <- out$sigma2/out$lambda
        out$rhohat <- NA
        out$shat.MLE <- NA
        out$best.model <- c(out$lambda, out$sigma2)
        out$warning <- "Maximum likelihood estimates not found with knots \n"
        print(paste(out$warning, " \n "))
    }
    if (!return.matrices) {
        out$xM <- NULL
        out$YM <- NULL
        out$x <- NULL
        out$y <- NULL
        out$matrices <- NULL
        out$weights <- NULL
    }
    if (!return.X) 
        out$matrices$X <- NULL
    out
}
