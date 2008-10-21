# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"Krig.engine.default" <-
function(out, verbose=FALSE){
# 
# matrix decompositions for computing estimate

#
# Computational outline:( "." is used for subscript)
#
# The form of the estimate is
#    fhat(x) = sum phi.j(x) d.j  + sum psi.k(x) c.k
#
# the {phi.j} are the fixed part of the model usually low order polynomials
# and is also referred to as spatial drift. 
#
# the {psi.k} are the covariance functions evaluated at the unique observation
# locations or "knots".  If xM.k is the kth unique location psi.k(x)= k(x, xM.k)
# xM is also out$knots in the code below.
#
# the goal is find decompositions that facilitate rapid solution for
# the vectors d and c. The eigen approach below was identified by
# Wahba, Bates Wendelberger and is stable even for near colinear covariance
# matrices. 
# This function does the main computations leading to the matrix decompositions.
# With these decompositions the coefficients of the solution are found in 
# Krig.coef and the GCV and REML functions in Krig.gcv.
#

#  First is an outline calculations with equal weights
#  T the fixed effects regression matrix  T.ij = phi.j(xM.i)
#  K the covariance matrix for the unique locations
# From the spline literature the solution solves the well known system 
# of two eqautions:

#    -K( yM - Td - Kc) + lambda *Kc = 0
#                 -T^t ( yM-Td -Kc) = 0
#
# Mulitple through by K inverse and substitute, these are equivalent to
#  
#  -1-   -( yM- Td - Kc) + lambda c = 0
#  -2-                        T^t c = 0       
#
#
#  A QR decomposition is done for   T= (Q.1,Q.2)R
#   by definition  Q.2^T T =0
#  
#  equation  -2- can be thought of as a constraint
# with  c= Q.2 beta2
# substitute in  -1-  and multiply through by Q.2^T
#
#      -Q.2^T yM  + Q.2^T K Q.2 beta2  + lambda beta2 = 0
#
#   Solving
#   beta2 = {Q.2^T K Q.2 + lambda I )^ {-1} Q.2^T yM
# 
# and so one sloves this linear system for beta2 and then uses
#     c= Q.2 beta2
#   to determine c.
# 
#  eigenvalues and eigenvectors are found for M= Q.2^T K Q.2
#     M = V diag(eta) V^T  
#  and these facilitate solving this system efficiently for
#  many different values of lambda.
#  create eigenvectors, D = (0, 1/eta)
#  and G= ( 0,0) %*% diag(D)
#         ( 0,V)
# so that
#
#          beta2 = G%*% ( 1/( 1+ lambda D)) %*% u
# with
#
#          u = (0, V Q.2^T W2 yM)
#
# Throughout keep in mind that M has smaller dimension than G due to
# handling the null space. 
#
# Now solve for d.
#  
# From -1-  Td = yM - Kc - lambda c
#      (Q.1^T) Td =  (Q.1^T) ( yM- Kc)    
#
#   ( lambda c is zero by -2-)
#   
#   so Rd = (Q.1^T) ( yM- Kc)
# use qr functions to solve triangular system in R to find d. 
#
#---------------------------------------------------------------------- 
# What about errors with a general precision matrix, W?
#  
# This is an important case because with replicated observations the
# problem will simplify into a smoothing problem with the replicate group
# means and unequal measurement error variances.
#
# the equations to solve are
#     -KW( yM - Td - Kc) + lambda *Kc = 0
#     -T^t W( yM-Td -Kc) =0
#
# Multiple through by K inverse and substitute, these are equivalent to 
#
#  -1b-      -W( yM- Td - Kc) + lambda c = 0
#  -2b-      (WT)^t c = 0       
#
# Let W2 be the symmetric square root of W,  W= W2%*% W2
# and W2.i be the inverse of W2.  
#
#  -1c-      -(  W2 yM - W2 T d - (W2 K W2) W2.ic) + lambda W2.i c = 0
#  -2c-      (W2T)^t  W2c = 0       
  
  
        Tmatrix<- do.call(out$null.function.name, 
                           c(out$null.args, list(x=out$xM,Z=out$ZM))  )
        if( verbose){
             cat(" Model Matrix: spatial drift and Z", fill=TRUE)
             print( Tmatrix)
        }

# Tmatrix premultiplied by sqrt of wieghts
        Tmatrix<- out$W2%d*%Tmatrix 

        qr.T <- qr(Tmatrix )

#
#verbose block
        if (verbose) {
            cat( "first 5 rows of qr.T$qr",fill=TRUE)
            print(qr.T$qr[1:5,])
        }
#
# find  Q_2 K Q_2^T  where K is the covariance matrix at the knot points
#
        
    tempM<-   t(
                 out$W2 %d*% 
                 do.call(out$cov.function.name, 
                 c(out$args, list(x1 = out$knots, x2 = out$knots))))
    tempM <- out$W2 %d*% tempM
    tempM <- qr.yq2(qr.T, tempM)
    tempM <- qr.q2ty(qr.T, tempM)

    np <- nrow(out$knots)
    nt <- (qr.T$rank)

    if (verbose) {
        cat("np, nt", np, nt, fill = TRUE)}

#
# Full set of decompositions for 
# estimator for nonzero lambda

            tempM <- eigen(tempM, symmetric=TRUE)
            D <- c(rep(0, nt), 1/tempM$values)
#
# verbose block
            if (verbose) {
                cat("eigen values:", fill = TRUE)
                print(D)
            }
#
# Find the transformed data vector used to 
# evaluate the solution, GCV, REML  at different lambdas
#
            u <- c( rep(0, nt),
                t(tempM$vectors) %*% 
                    qr.q2ty(qr.T, c(out$W2%d*%out$yM ) ))
#
#
            return( list(D = D,  qr.T = qr.T, 
                decomp = "WBW", V = tempM$vectors,u=u, nt=nt, np=np))
}

