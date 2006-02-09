"Krig" <-
function (x, Y, cov.function = "exp.cov", lambda = NA, df = NA, 
    cost = 1, knots, weights = NULL, m = 2, 
return.matrices = TRUE, 
    nstep.cv = 80, scale.type = "user", x.center = rep(0, ncol(x)), 
    x.scale = rep(1, ncol(x)), rho = NA, sigma2 = NA, method = "GCV", 
    decomp = "DR", verbose = FALSE, cond.number = 10^8, mean.obj = NULL, 
    sd.obj = NULL, yname = NULL, return.X = TRUE, null.function = make.tmatrix, 
    offset = 0, outputcall = NULL, cov.args = NULL, na.rm = FALSE, 
    ...) 
{

# setup some components of Krig output object
    out <- list()
    class(out) <- c("Krig")
    if (!is.character(cov.function)) {
        if (is.function(cov.function)) 
            cov.function <- as.character(substitute(cov.function))
    }
#    name of covariance function
    out$call.name <- cov.function

    if (is.null(outputcall)) {
        out$call <- match.call()
    }
    else {
        out$call <- outputcall
    }

# check for missing values in y

    if (sum(is.na(Y)) != 0 & !na.rm) {
        stop("Need to remove missing values or use: na.rm=TRUE in the call")
    }
# define weights
# name of Y data

    if (is.null(yname)) {
        out$yname <- deparse(substitute(Y))
    }
    else {
        out$yname <- yname
    }

# various things to save 
    out$decomp <- decomp
    out$make.tmatrix <- null.function
    out$offset <- offset
    out$cost <- cost

# save extra arguments to the covariance function 

    if (!is.null(cov.args)) 
        out$args <- cov.args
    else out$args <- list(...)
    if (verbose) {
        cat(" Local cov function arguments ", fill = TRUE)
        print(out$args)
        cat(" covariance function used is : ", fill = TRUE)
        print(out$call.name)
    }

# coerce x to be a matrix and figure out the column names

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

# coerce Y to be a vector 


    Y <- c(Y)

# setup default for weights
    if( is.null( weights)) weights<- rep( 1,length(Y))

# excise the NA's 
    if (na.rm) {
        ind <- is.na(Y)
        if (any(ind)) {
            Y <- Y[!ind]
            x <- x[!ind, ]
            weights<- weights[!ind]
            warning("NA's have been removed from Y 
                and the corresponding X rows .")
        }
    }
  if ( verbose){
   print( Y)
   print( x)
   print( weights)
  }
# save basic info 
    N <- length(Y)
    out$N <- N
    out$y <- Y
    out$x <- x
    out$weights <- weights

# standardized Y if this is a correlation model
    if (!is.null(sd.obj) & !is.null(mean.obj)) {
        correlation.model <- TRUE
    }
    else correlation.model <- FALSE
    out$correlation.model <- correlation.model
    if (correlation.model) {
        Yraw <- Y
        out$mean.obj <- mean.obj
        out$sd.obj <- sd.obj
        Y <- (Y - predict(mean.obj, x))/predict(sd.obj, x)
        if (verbose)
            print(N) 
            print(Y)
    }

# figure out if there are replicates and if so collapse onto 
# means, unique locations and average weights. 

    rep.info <- cat.matrix(x)
    out$rep.info <- rep.info

    if (verbose) {
        cat("rep.info", fill = TRUE)
        print(rep.info)
    }
    if (max(rep.info) == N) {
        shat.rep <- NA
        shat.pure.error <- NA
        out$pure.ss <- 0
        YM <- Y
        weightsM <- weights
        xM <- as.matrix(x[!duplicated(rep.info), ])
    }
    else {
        rep.info.aov <- fast.1way(rep.info, Y, weights)
        shat.pure.error <- sqrt(rep.info.aov$MSE)
        shat.rep <- shat.pure.error
        YM <- rep.info.aov$means
 
        weightsM <- rep.info.aov$w.means
        xM <- as.matrix(x[!duplicated(rep.info), ])
        out$pure.ss <- rep.info.aov$SSE

        if (verbose) {
            cat(" rep info", fill = TRUE)
            print(rep.info.aov)
            cat( "weightsM", fill=TRUE)
            print( weightsM)
        }
    }
# save replicate info to output object 
   out$yM <- YM
    out$xM <- xM
    out$weightsM <- weightsM
    out$shat.rep <- shat.rep
    out$shat.pure.error <- shat.pure.error

# default choice for knots are the unique locations
# set flags whether an MLE for lambda is possible

    if (missing(knots)) {
        knots <- xM
        mle.calc <- TRUE
        knot.model <- FALSE
    }
    else {
        mle.calc <- FALSE
        knot.model <- TRUE
    }

    out$knot.model <- knot.model
    out$mle.calc <- mle.calc
    knots <- as.matrix(knots)
    out$knots <- knots
#
# scale locations transform both x and the knots

    xM <- transformx(xM, scale.type, x.center, x.scale)
    transform <- attributes(xM)

    if (verbose) {
        cat("transform", fill = TRUE)
        print(transform)
    }

    knots <- scale(knots, center = transform$x.center, scale = transform$x.scale)

    if (verbose) {
        cat("knots in transformed scale", fill = TRUE)
        print(knots)
    }
    out$transform <- transform

#
#  setup method for finding lambda, rho and sigma

    if (!is.na(lambda) | !is.na(df)) 
        method <- "user"
    if (!is.na(rho) & !is.na(sigma2)) {
        lambda <- sigma2/rho
        method <- "user"
    }
    just.solve <- (lambda[1] == 0)
    if (is.na(just.solve)) 
        just.solve <- FALSE
    if (verbose) 
        cat("lambda", lambda, fill = TRUE)

# set dimension of locations
    d <- ncol(xM)

# begin main computation block there are two main branches
#  either the DR method or WBW method
   
 if (decomp == "DR") {
        qr.T <- qr(out$make.tmatrix(knots, m))
#        if (verbose) {
#            print(qr.T)
#        }
        tempM <- qr.yq2(qr.T, do.call(out$call.name, c(out$args, 
            list(x1 = knots, x2 = knots))))
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
        tempM <- sqrt(weightsM) * t(sqrt(weightsM) * t(do.call(out$call.name, 
            c(out$args, list(x1 = knots, x2 = knots)))))
        tempM <- qr.yq2(qr.T, tempM)
        tempM <- qr.q2ty(qr.T, tempM)
    }
    np <- nrow(knots)
    nt <- (qr.T$rank)
    out$np <- np
    out$nt <- nt
    if (verbose) 
        cat("np, nt", np, nt, fill = TRUE)
    if (verbose) 
        print(knots)
    if (just.solve) {
        beta <- qr.coef(qr(cbind(out$make.tmatrix(xM, m), qr.yq2(qr.T, 
            do.call(out$call.name, c(out$args, list(x1 = xM, 
                x2 = knots)))))), YM)
    }
    else {
        if (decomp == "DR") {
            if (verbose) 
                cat("Type of decomposition", decomp, fill = TRUE)
            H <- matrix(0, ncol = np, nrow = np)
            H[(nt + 1):np, (nt + 1):np] <- tempM
            X <- cbind(out$make.tmatrix(xM, m), qr.yq2(qr.T, 
                do.call(out$call.name, c(out$args, list(x1 = xM, 
                  x2 = knots)))))
            if (verbose) {
                cat( "weightsM", fill=TRUE)
                print(weightsM)
                cat("first 3 rows of X", fill = TRUE)
                #print(X[1:3, ])
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
                cat("singular values:", fill = TRUE)
                print(D)
            }
            D[(1:nt) + (np - nt)] <- 0
            G <- B %*% U
            u <- t(G) %*% t(X) %*% (weightsM * YM)
            if (verbose) {
                cat("DR: u", fill = TRUE)
                print(u)
            }
            out$pure.ss <- sum(weightsM * (YM - X %*% G %*% u)^2) + 
                out$pure.ss
            if (verbose) {
                cat("pure.ss", fill = TRUE)
                print(out$pure.ss)
            }
            out$matrices <- list(B = B, U = U, u = u, D = D, 
                G = G, qr.T = qr.T, X = X)
        }
        if (decomp == "WBW") {
            temp <- svd(tempM)[c("d", "v")]
            D <- c(rep(0, nt), 1/temp$d)
            if (verbose) {
                cat("singular values:", fill = TRUE)
                print(D)
            }
            G <- matrix(0, ncol = np, nrow = np)
            G[(nt + 1):np, (nt + 1):np] <- temp$v
            G <- G * matrix(D, ncol = np, nrow = np, byrow = TRUE)
            u <- c(rep(0, nt), t(temp$v) %*% qr.q2ty(qr.T, sqrt(weightsM) * 
                YM))
            if (verbose) {
                cat("WBW: u", fill = TRUE)
                print(u)
            }
            if (verbose) {
                cat("WBW: pure.ss", fill = TRUE)
                print(out$pure.ss)
            }
            out$matrices <- list(u = u, D = D, G = G, qr.T = qr.T, 
                decomp = decomp, V = temp$v)
        }
        if (verbose) {
            cat("call to gcv minimization", fill = TRUE)
        }
        gcv.out <- gcv.Krig(out, nstep.cv = nstep.cv, verbose = verbose, 
            cost = out$cost, offset = out$offset, lambda = lambda)
        gcv.grid <- gcv.out$gcv.grid
        out$gcv.grid <- gcv.grid
        out$lambda.est <- gcv.out$lambda.est
        if (verbose) {
            print(out$gcv.grid)
        }
    }
    out$method <- method
    if (method == "user") {
        if (!is.na(df)) 
            lambda <- Krig.df.to.lambda(df, D)
        lambda.best <- lambda
        out$lambda.est.user <- summary.gcv.Krig(out, lambda.best, 
            offset = out$offset, cost = out$cost)
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
            cat("trace of A", fill = TRUE)
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
        cat(names(out), fill = TRUE)
    }
    out$fitted.values <- predict.Krig(out, x, eval.correlation.model = FALSE)
    out$residuals <- Y - out$fitted.values
    if (verbose) {
        cat("resdiuals", out$residuals, fill = TRUE)
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

