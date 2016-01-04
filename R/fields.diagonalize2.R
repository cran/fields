# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"fields.diagonalize2" <- function(A, B, verbose = FALSE) {
    M <- nrow(A)
    hold.AB <- eigen(A + B, symmetric = TRUE)
    if (verbose) {
        cat("log 10 condition number of A +B  in fields.diagonlize2", 
            fill = TRUE)
        print(log10(max(hold.AB$values)/min(hold.AB$values)))
    }
    #   inverse square root of  A+B
    hold.AB <- (t(hold.AB$vectors) * (1/sqrt(hold.AB$values)))
    hold.B <- eigen(hold.AB %*% A %*% t(hold.AB), symmetric = TRUE)
    G <- t(hold.B$vectors) %*% hold.AB
    D.A <- hold.B$values
    # remove some large temporary matrices.
    remove(hold.AB)
    remove(hold.B)
    # crank on finding G and D.
    G <- (1/sqrt(D.A)) * G
    D <- colSums(t(G) * (B) %*% t(G))
    # sort from largest to smallest  and take transpose ---
    #  this will now matches old version in fields.diagonalize
    D <- D[M:1]
    G <- t(G[M:1, ])
    # to test:
    #    test.for.zero(  t(G) %*% (A) %*% (G), diag(1,M), tag='A test' )
    #    test.for.zero(  t(G) %*% (B) %*% (G), diag(D,M), tag='B test' )
    list(G = G, D = D)
}
