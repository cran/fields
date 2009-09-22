# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
spind2full <- function(obj) {
    # create empty matrix and stuff at non zero locations
    temp <- matrix(0, obj$da[1], obj$da[2])
    temp[obj$ind] <- obj$ra
    return(temp)
}
spind2spam <- function(obj) {
    # sort on rows and then columns to make sure they are in order
    ii <- order(obj$ind[, 1], obj$ind[, 2])
    #    shuffle indices and entries so they are in row order
    obj$ind <- obj$ind[ii, ]
    obj$ra <- obj$ra[ii]
    ia <- obj$ind[, 1]
    # define total number of nonzero elements
    M <- length(ia)
    # find places where rows change
    hold <- diff(c(0, ia, M + 1))
    # Note: 1:M is the cumsum for elements.
    ia <- (1:(M + 1))[hold != 0]
    # check if there is a missing row.
    # if so stop -- because bad things happen ..
    if (length(unique(obj$ind[, 1])) < obj$da[1]) {
        # The offending rows -- concise but inpentrable R code!
        ind.missing <- (1:obj$da[1])[-unique(obj$ind[, 1])]
        stop(paste("Row(s)", ind.missing, "  missing in matrix"))
    }
    return(new("spam", entries = obj$ra, colindices = obj$ind[, 
        2], rowpointers = ia, dimension = obj$da))
}
spam2spind <- function(obj) {
    # diff gives the number of nonzero elements in each row
    I <- rep((1:obj@dimension[1]), diff(obj@rowpointers))
    list(ind = cbind(I, obj@colindices), da = obj@dimension, 
        ra = obj@entries)
}
spam2full <- function(obj) {
    spind2full(spam2spind(obj))
}
