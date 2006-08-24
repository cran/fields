"Krig.coef" <-
function (out, lambda = out$lambda, y = NULL, yM=NULL, verbose=FALSE) 
{

#
# NOTE default value of lambda used from Krig object.

#
# Determine whether to collapse onto means of replicates ( using y)
# or if data as replicate means that have been passed (as yM) use it!
# If both y and YM are null then just use out$yM 
# for readbality of this function all this torturous logic happens in 
#  Krig.ynew. 

    out2 <- Krig.ynew(out, y, yM)
    temp.yM <- out2$yM

    nt <- out$nt
    np <- out$np
    u <- NA
    call.name<- out$cov.function.name

    if (out$decomp == "DR") {
      
    # X is the monster matrix ...

      X <- cbind(
        out$make.tmatrix(out$xM, out$m),
        do.call(call.name, c(out$args, list(x1 = out$xM, x2 =out$knots)))
           )

        u <- t(out$matrices$G) %*% t(X) %*% (out$weightsM * 
            temp.yM)
        beta <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$D)) * 
            u)
        temp.d <- beta[1:nt]
        temp.c <- beta[(nt + 1):np]
        temp <- X %*% out$matrices$G %*% u
        temp <- sum(out$weightsM * (temp.yM - temp)^2)
        out2$pure.ss <- temp + out2$pure.ss
    }
    if (out$decomp == "WBW") {
        u <- c(rep(0, out$nt), t(out$matrices$V) %*% qr.q2ty(out$matrices$qr.T, 
            sqrt(out$weightsM) * temp.yM))
        beta <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$D)) * 
            u)
        temp.c <- c(qr.qy(out$matrices$qr.T, c(rep(0, nt), beta[(nt + 
            1):np])))
        temp.c <- temp.c * sqrt(out$weightsM)
        temp <- temp.yM - lambda * temp.c - do.call(call.name, 
            c(out$args, list(x1 = out$knots, x2 = out$knots, 
                C = temp.c)))
        temp <- sqrt(out$weightsM) * temp
        temp.d <- qr.coef(out$matrices$qr.T, temp)
    }
    if (out$decomp == "cholesky") {
        if (lambda != out$matrices$lambda) {
            stop("New lambda can not be used with cholesky decomposition")
        }
        Tmatrix <- out$make.tmatrix(out$knots, out$m)
        temp.d <- qr.coef(out$matrices$qr.VT, 
         forwardsolve(out$matrices$Mc, 
            transpose = TRUE, temp.yM, upper.tri=TRUE))
        temp.c <- forwardsolve(out$matrices$Mc, transpose = TRUE, 
            temp.yM - Tmatrix %*% temp.d, upper.tri=TRUE)
        temp.c <- backsolve(out$matrices$Mc, temp.c)
    }
    if (out$decomp == "cholesky.knots") {
        if (lambda != out$matrices$lambda) {
            stop("New lambda can not be used with cholesky decomposition")
        }
    # form K matrix
         K<-  do.call( call.name,
            c(out$args, list(x1 = out$xM, x2 = out$knots) ) )
#

        Tmatrix <- out$make.tmatrix(out$xM, out$m)

        wY <- out$weightsM * temp.yM
        temp0 <- t(K) %*% (out$weightsM * Tmatrix)
        temp1 <- forwardsolve(out$matrices$Mc, temp0, transpose = TRUE,
                    upper.tri=TRUE)
        qr.Treg <- qr(t(Tmatrix) %*% (out$weightsM * Tmatrix) - 
            t(temp1) %*% temp1)
        temp0 <- t(K) %*% wY
        temp3 <- t(Tmatrix) %*% wY - t(temp1) %*% forwardsolve(out$matrices$Mc, 
            temp0, transpose = TRUE, upper.tri=TRUE)
        temp.d <- qr.coef(qr.Treg, temp3)
        temp1 <- t(K) %*% (wY - out$weightsM * (Tmatrix) %*% 
            temp.d)
        temp.c <- forwardsolve(out$matrices$Mc, transpose = TRUE, 
            temp1, upper.tri=TRUE)
        temp.c <- backsolve(out$matrices$Mc, temp.c)
    }
    return(c(out2, list(c = c(temp.c), d = c(temp.d), u = c(u))))
}

