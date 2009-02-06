# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"vgram" <-
function (loc, y, id = NULL, d = NULL, lon.lat = FALSE, dmax = NULL, 
    N = NULL, breaks = NULL) 
{
# coerce to matrix
    y<- cbind( y)

# if nearest neighbor indices are missing create all possible pairs. 
    if (is.null(id)) {
        n <- nrow(loc)
        ind <- rep(1:n, n) > rep(1:n, rep(n, n))
        id <- cbind(rep(1:n, n), rep(1:n, rep(n, n)))[ind, ]
    }

# if distances are missing calculate these

    if (is.null(d)) {
        loc <- as.matrix(loc)
        if (lon.lat) {
            d <- rdist.earth(loc)[id]}
        else {
            d <- rdist(loc, loc)[id]}
    }

#    
# calculating variogram will average over columns if y is a matrix
# 
    vg <- 0.5 * rowMeans( cbind((y[id[, 1],] - y[id[, 2],])^2) , na.rm=TRUE)

#
#information for returned object
#
    call <- match.call()
     if (is.null(dmax)) {
        dmax <- max(d)}

    od <- order(d)
    d <- d[od]
    vg <- vg[od]
    ind <- d <= dmax & !is.na(vg)

## add a binned  variogram if breaks are supplied
    out <- list(d = d[ind], vgram = vg[ind], call = call)

    if (!is.null(breaks) | !is.null(N)) {
        out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks)) }

    out
}

