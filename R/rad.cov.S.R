"rad.cov.S" <-
function (x1, x2, p = 1, with.log = TRUE, with.constant = TRUE) 
{
    if (!is.matrix(x1)) 
        x1 <- as.matrix(x1)
    if (!is.matrix(x2)) 
        x2 <- as.matrix(x2)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    m <- (d + p)/2
    temp <- (outer(x1[, 1], x2[, 1], "-"))^2
    if (ncol(x1) > 1) {
        for (k in 2:ncol(x1)) {
            temp <- temp + (outer(x1[, k], x2[, k], "-"))^2
        }
    }
    if ((d%%2 == 0) & (with.log)) 
        temp <- ifelse(temp < 1e-10, 0, temp^(p/2) * log(temp))
    else temp <- temp^(p/2)
    if (with.constant) {
        Amd <- radbas.constant(m, d)
    }
    else {
        Amd <- 1
    }
    Amd * matrix(temp, ncol = n2, nrow = n1)
}
