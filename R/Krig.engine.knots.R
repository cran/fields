# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"Krig.engine.knots" <- function(out, verbose = FALSE) {
    #
    # matrix decompostions for computing estimate when
    # knots are present
    # QR decomposition of null space regression matrix
    Tmatrix <- do.call(out$null.function.name, c(out$null.args, 
        list(x = out$xM, Z = out$ZM)))
    qr.T <- qr(c(sqrt(out$weightsM)) * Tmatrix)
    nt <- ncol(Tmatrix)
    np <- nrow(out$knots) + nt
    if (verbose) {
        cat(nt, np, fill = TRUE)
    }
    # H is the penalty matrix in the ridge regression format
    # first part is zero because no penalty on part of estimator
    # spanned by T matrix
    H <- matrix(0, ncol = np, nrow = np)
    H[(nt + 1):np, (nt + 1):np] <- do.call(out$cov.function.name, 
        c(out$args, list(x1 = out$knots, x2 = out$knots)))
    # X is the monster ...
    X <- cbind(do.call(out$null.function.name, c(out$null.args, 
        list(x = out$xM, Z = out$ZM))), do.call(out$cov.function.name, 
        c(out$args, list(x1 = out$xM, x2 = out$knots))))
    if (verbose) {
        cat("first lines of X", fill = TRUE)
        print(X[1:5, ])
    }
    #    sqrt(weightsM) * X
    XTwX <- t(X * out$weightsM) %*% X
    #
  # then  B= G(I-D)G^T
  # New version of diagonalize may be more stable 
    out2 <- fields.diagonalize2((XTwX), H)
    D <- out2$D
    if (verbose) {
        cat("D;", fill = TRUE)
        cat(out2$D, fill = TRUE)
    }
#
#  G should satisfy:
#     t(G) %*% XTwX %*%G = I and  t(G)%*%H%*%G = D
#
#     and 
#      solve( XtwX + lambda H) =  G%*%diag( 1/(1+ lambda*D))%*%t(G)
# 

#  save XG to avoid an extra multiplication.
    XG<- X%*% out2$G 

    u <- t(XG) %*% (out$weightsM * out$yM)
    #
    # adjust pure sum of squares to be that due to replicates
    # plus that due to fitting all the basis functions without
    # any smoothing. This will be the part of the RSS that does not
    # change as lambda is varied ( see e.g. gcv.Krig)
    #
    pure.ss <- sum(out$weightsM * (out$yM - XG %*% 
        u)^2) + out$pure.ss
    if (verbose) {
        cat("total pure.ss from reps, reps + knots ", fill = TRUE)
        print(out$pure.ss)
        print(pure.ss)
    }

    #
    # in this form  the solution is (d,c)= G( I + lambda D)^-1 u
    # fitted.values = X ( d,c)
    #
    # output list
# last D eigenvalues are zero due to null space of penalty
# OLD code:    D[(np - nt + 1):np] <- 0
# this should be enforced to machine precision from diagonalization. 


    list(u = u, D = D, G = out2$G, qr.T = qr.T, decomp = "DR", 
        nt = nt, np = np, pure.ss = pure.ss)
}
