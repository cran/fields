"predict.tps" <-
function (out, x, y, lambda, df, omega, u = NA, derivative = 0, 
    model = NA) 
{
    if (!is.na(model)) 
        lambda <- model
    if (out$tag != 1) 
        stop("This is an old tps object please rerun\n\ntps to get right coefficients")
    if (missing(x)) 
        x <- out$x
    x <- as.matrix(x)
    xc <- out$transform$x.center
    xs <- out$transform$x.scale
    n <- nrow(x)
    p <- ncol(x)
    knots <- scale(out$knots, xc, xs)
    x <- scale(x, xc, xs)
    nt <- out$nt
    np <- out$np
    dtemp <- out$d
    ctemp <- out$c
    if (!missing(lambda) | !missing(df) | !missing(y)) {
        if (missing(lambda)) 
            lambda <- out$lambda
        if (!missing(df)) 
            lambda <- tps.df.to.lambda(df, out$matrices$D)
        if (!missing(y)) {
            if (is.null(out$matrices$X)) {
                stop("big X matrix not part of matrices\nobject from tps!")
            }
            u <- t(out$matrices$X %*% out$matrices$G) %*% (y * 
                out$weights)
        }
        else {
            u <- out$matrices$u
        }
        omega <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$D)) * 
            u)
        dtemp <- omega[1:nt]
        temp <- c(rep(0, nt), omega[(nt + 1):np])
        ctemp <- c(qr.qy(out$matrices$qr.T, temp))
    }
    if (!missing(omega)) {
        dtemp <- omega[1:nt]
        temp <- c(rep(0, nt), omega[(nt + 1):np])
        ctemp <- c(qr.qy(out$matrices$qr.T, temp))
    }
    if (derivative == 0) {
        return(make.tmatrix(x, out$m) %*% dtemp + make.Kc(x, 
            knots, ctemp, p = out$power, with.constant = out$with.constant))
    }
    if (derivative == 1) {
        temp <- matrix(1/xs, ncol = p, nrow = n, byrow = T)
        return((make.DTd(x, dtemp, m = out$m) + make.DKc(x, knots, 
            ctemp, p = out$power, with.constant = out$with.constant)) * 
            temp)
    }
}
