"rdist.S" <-
function (x1, x2, lon.lat = FALSE) 
{
    if (!is.matrix(x1)) {
        x1 <- as.matrix(x1)
    }
    if (!is.matrix(x2)) {
        x2 <- as.matrix(x2)
    }
    if (lon.lat) {
        rdist.earth(x1, x2)
    }
    else {
        temp <- (outer(x1[, 1], x2[, 1], "-"))^2
        if (ncol(x1) > 1) {
            for (k in 2:ncol(x1)) {
                temp <- temp + (outer(x1[, k], x2[, k], "-"))^2
            }
        }
        sqrt(temp)
    }
}
