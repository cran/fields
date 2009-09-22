# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"plot.vgram.matrix" <- function(x, ...) {
    ind <- x$ind
    ir <- range(ind[, 1])
    jr <- range(ind[, 2])
    # x and y grid values
    temp.list <- list(x = (ir[1]:ir[2]) * x$dx, y = (jr[1]:jr[2]) * 
        x$dy)
    # fill in a matrix with variogram values
    ind2 <- cbind(ind[, 1] - min(ind[, 1]) + 1, ind[, 2] - min(ind[, 
        2]) + 1)
    temp <- matrix(NA, nrow = max(ind2[, 1]), ncol = max(ind2[, 
        2]))
    temp[ind2] <- x$vgram.full
    temp.list$z <- temp
    # plot it!
    image.plot(temp.list, ...)
}
