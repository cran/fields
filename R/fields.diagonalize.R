# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"fields.diagonalize" <- function(A, B) {

    hold <- eigen(A, symmetric = TRUE)
    # square root of A
    hold2 <- (t(hold$vectors) * sqrt(1/hold$values))
    #
    # A.inv.sqrt = hold2
    # A.inv = hold%*% t(hold2)
    #
    # eigen decomp of  A.inv.sqrt B t( A.inv.sqrt)
    #
    hold <- eigen((hold2) %*% B %*% t(hold2), symmetric = TRUE)
    # the magic G matrix used throughout fields.
    G <- t(hold2) %*% hold$vectors
    #
    # Note:
    # G simultaneously diagonalizes two matrices:
    #
    # G^T A G= I
    # G^T B G= D
    #
    # and in terms of application we also have the useful
    # diagonalization
    #
    #  (A +lambda B)^{-1} =  G( I + lambda D)^{-1} G^T
    list(G = G, D = hold$values)
}
