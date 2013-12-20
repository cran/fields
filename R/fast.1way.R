# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"fast.1way" <- function(lev, y, w = rep(1, length(y))) {
    # w are proportional to reciprocal variance.
    if (!is.matrix(y)) {
        y <- as.matrix(y)
    }
    N <- nrow(y)
    NC <- ncol(y)
    # ordered unique values of lev
    tags <- lev[!duplicated(lev)]
    NR <- length(tags)
    # lev are now integer tags
    lev <- match(lev, tags)
    #
    means <- matrix(NA, nrow = NR, ncol = NC)
    # add together weights with same lev
    w.means <- c(tapply(w, lev, sum))
    for (k in 1:NC) {
        # find weighted means for each lev
        means[, k] <- (tapply(y[, k] * w, lev, sum)/w.means)
    }
    # find SS
    SSE <- colSums((w * (y - means[lev, ])^2))
    MSE <- SSE/(N - NR)
    list(n = N, means = means, SSE = SSE, w.means = w.means, 
        MSE = MSE, lev = lev, tags = tags)
}
