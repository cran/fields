"make.tmatrix.S" <-
function (x, m = 2) 
{
    if (!is.matrix(x)) 
        x <- as.matrix(x)
    p <- m - 1
    d <- ncol(x)
    n <- nrow(x)
    ptab <- NA
    ind <- (1:d) + 1
    end <- ind[d]
    N <- as.integer(prod(((p + 1):(p + d))/(1:d)))
    A <- matrix(NA, ncol = N, nrow = n)
    ptab <- matrix(0, ncol = N, nrow = d)
    unit <- diag(rep(1, d))
    A[, 1] <- rep(1, n)
    if (p == 0) {
        attr(A, "ptab") <- t(ptab)
        return(A)
    }
    A[, ind] <- x
    ptab[, ind] <- unit
    if (p > 1) {
        for (j in 2:p) {
            k1 <- end + 1
            for (i in 1:d) {
                k2 <- end - ind[i] + k1
                A[, k1:k2] <- x[, i] * A[, ind[i]:end]
                ptab[, k1:k2] <- ptab[, ind[i]:end] + unit[, 
                  i]
                ind[i] <- k1
                k1 <- k2 + 1
            }
            end <- k2
        }
    }
    attr(A, "ptab") <- t(ptab)
    A
}
