"qsreg" <-
function (x, y, lam = NA, maxit = 50, maxit.cv = 10, tol = 1e-07, 
    offset = 0, sc = sqrt(var(y)) * 1e-05, alpha = 0.5, wt = rep(1, 
        length(x)), cost = 1, nstep.cv = 100, hmin = NA, hmax = NA, 
    trmin = 2 * 1.05, trmax = 0.95 * length(unique(x))) 
{
    if (!is.loaded(symbol.For("css"))) {
        temp <- dyn.load(paste(FIELDS.BIN, "fields.o", sep = ""), 
            2)
    }
    out <- list()
    class(out) <- c("qsreg")
    N <- length(y)
    out$N <- N
    xgrid <- sort(unique(x))
    if (length(x) != length(y)) 
        stop(" X and Y do not match")
    if (is.na(lam[1])) 
        do.cv <- T
    else do.cv <- F
    if (do.cv) {
        if (is.na(hmin)) {
            hmin <- 0
            for (k in 1:25) {
                b <- qsreg.trace(lam = as.double(exp(hmin)), 
                  x = x, y = y, wt = wt, cost = cost, maxit = maxit, 
                  tol = tol, sc = sc, alpha = alpha)
                if (b > trmax) {
                  break
                }
                hmin <- hmin - 1
            }
        }
        if (is.na(hmax)) {
            hmax <- 0
            for (k in 1:25) {
                b <- qsreg.trace(lam = as.double(exp(hmax)), 
                  x = x, y = y, wt = wt, cost = cost, maxit = maxit, 
                  tol = tol, sc = sc, alpha = alpha)
                if (b < trmin) {
                  break
                }
                hmax <- hmax + 1
            }
        }
        a <- .Fortran("cvrcss", n = as.integer(N), x = as.double(x), 
            y = as.double(y), wt = as.double(wt), sy = as.double(rep(0, 
                N)), diag = as.double(rep(0, N)), din = as.double(c(cost, 
                offset, maxit, tol, sc, alpha)), dout = as.double(rep(0, 
                4)), nstep = as.integer(nstep.cv), maxit = as.integer(maxit.cv), 
            trmin = as.double(trmin), trmax = as.double(trmax), 
            hmin = as.double(hmin), hmax = as.double(hmax), hopt = as.double(-1), 
            vopt = as.double(-1), tropt = as.double(-1), mxstep = as.integer(nstep.cv), 
            tabout = as.double(rep(0, 4 * nstep.cv)), ierr = as.integer(0))
        if (a$ierr == -1) {
            cat("minimum CV is at the\nboundary of the grid\nfor minimization", 
                fill = T)
        }
        if (a$ierr > 0) {
warning("Plot the returned object to examine the CV function. Refined CV  search failed  ")
        }
    }
    if (do.cv) {
        lam <- exp(a$hopt)
        wt <- a$wt
    }
    b <- list()
    NL <- length(lam)
    NG <- length(xgrid)
    h <- log(lam)
    residuals <- matrix(0, ncol = NL, nrow = N)
    diagA <- residuals
    cv <- rep(0, NL)
    predicted <- matrix(0, ncol = NL, nrow = NG)
    trace <- rep(0, NL)
    converge <- rep(0, NL)
    for (k in 1:NL) {
        b <- .Fortran("rcss", h = as.double(h[k]), npoint = as.integer(N), 
            x = as.double(x), y = as.double(y), wt = as.double(wt), 
            sy = as.double(rep(0, N)), trace = as.double(0), 
            diag = as.double(rep(0, N)), cv = as.double(0), ngrid = as.integer(NG), 
            xg = as.double(xgrid), yg = as.double(rep(0, NG)), 
            job = as.integer(c(3, 3, 0)), ideriv = as.integer(0), 
            din = as.double(c(cost, offset, maxit, tol, sc, alpha)), 
            dout = as.double(rep(0, 4)), ierr = as.integer(0))
        residuals[, k] <- y - b$sy
        diagA[, k] <- b$diag
        cv[k] <- b$dout[4]
        trace[k] <- b$trace
        predicted[, k] <- b$yg
        converge[k] <- b$dout[1]
    }
    if (do.cv) {
        cv.grid <- matrix(a$tabout, ncol = 4)
        cv.grid[, 1] <- exp(cv.grid[, 1])
    }
    else cv.grid <- cbind(lam, trace, cv, converge)
    dimnames(cv.grid) <- list(NULL, c("lambda", "trace", "CV", 
        "iterations"))
    out$call <- match.call()
    out$x <- x
    out$y <- y
    out$residuals <- residuals
    out$fitted.values <- y - residuals
    out$predicted <- list(x = xgrid, y = predicted)
    out$cv.grid <- cv.grid
    out$trace <- trace
    out$lambda <- lam
    out$diagA <- diagA
    out$sc <- sc
    out$alpha <- alpha
    out
}
