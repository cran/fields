"Krig.engine.default" <-
function(out, verbose=FALSE){
# 
# matrix decompositions for computing estimate

#
# Computational outline ( "." is used for subscript)
#
# The form of the estimate is
#    fhat(x) = sum phi.j(x) d.j  + sum psi.k c.k
#
# the {phi.j} are the fixed part of the model usually low order polynomials
# and is also referred to as spatial drift.
#
# the {psi.k} are the covariance function evaluated at the unique observation
# locations.  if xM.k is the kth unique location psi.k(x)= k(x, xM.k)
# xM is also out$knots in the code below.
#
# the goal is find decompositions that facilitate rapid solution for
# the vectors d and c. 
#

#  W  be the weight matrix for unique locations e.g. diag( out$weightsM)
#  T the fixed effects regression matrix  T.ij = phi.j(xM.i)
#  K the covariance matrix for the unique locations

# From the spline literature the solution solves the well known system 
# of two eqautions:

#    -KW( yM - Td - Kc) + lambda *Kc = 0
#     -T^t W ( yM-Td -Kc) =0

# Divide through by K and substitute, these are equivalent to 

#       -W( yM- Td - Kc) + lambda c = 0
#        T^t c = 0       



#
#  A QR decomposition is done for   sqrt(W)T= (Q.1,Q.2)R
#
#  eigenvalues and eigenvectors found for M= Q.2^T K Q.2 
#  


#
# QR decomposition for the matrix defining the fixed effects or the 
# null space. 
# Note premultiplication by sqrt of weights 
#   QR   where Q= ( Q_1, Q_2)  Q_1 spans column space of T

        qr.T <- qr(c(sqrt(out$weightsM)) * out$make.tmatrix(out$xM, out$m))
#
#verbose block
        if (verbose) {
cat( "first 5 rows of qr.T$qr",fill=TRUE)
            print(qr.T$qr[1:5,])
        }
#
# find  Q_2 K Q_2^T  where K is the covariance matrix at the knot points
#
        tempM <- sqrt(out$weightsM) * t(sqrt(out$weightsM) * 
            t(do.call(out$cov.function.name, 
            c(out$args, list(x1 = out$knots, x2 = out$knots)))))
        tempM <- qr.yq2(qr.T, tempM)
        tempM <- qr.q2ty(qr.T, tempM)

    np <- nrow(out$knots)
    nt <- (qr.T$rank)

    if (verbose) {
        cat("np, nt", np, nt, fill = TRUE)
}

#
# Full set of decompositions for 
# estimator for nonzero lambda

            temp <- eigen(tempM, symmetric=TRUE)
            D <- c(rep(0, nt), 1/temp$values)
#
# verbose block
            if (verbose) {
                cat("eigen values:", fill = TRUE)
                print(D)
            }
#
# Form the matrix decompositions and transformed data vector used to 
# evaluate the solution, GCV, REML  at different lambdas
#
            G <- matrix(0, ncol = np, nrow = np)
            G[(nt + 1):np, (nt + 1):np] <- temp$vectors
            G <- G * matrix(D, ncol = np, nrow = np, byrow = TRUE)
            u <- c(rep(0, nt), t(temp$vectors) %*% qr.q2ty(qr.T, 
                        sqrt(out$weightsM) * out$yM))
#
# verbose block
            if (verbose) {
                print(u)
                print(out$pure.ss)
            }
#
            return( list(u = u, D = D, G = G, qr.T = qr.T, 
                decomp = "WBW", V = temp$vectors, nt=nt, np=np))
}

