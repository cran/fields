# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
spind2full <- function(obj) {
    # create empty matrix and stuff at non zero locations
    temp <- matrix(0, obj$da[1], obj$da[2])
    temp[obj$ind] <- obj$ra
    return(temp)
}
spind2spam <- function(obj, add.zero.rows = TRUE) {
    # Check if there is one or more  missing rows. If so either stop or fill in these rows
    #   with a zero value
    rows.present <- unique(obj$ind[, 1])
    if (length(rows.present) < obj$da[1]) {
        # The missing row indices
        ind.missing <- (1:obj$da[1])[-rows.present]
        N.missing <- length(ind.missing)
        if (!add.zero.rows) {
            cat(N.missing, " missing row(s)", fill = TRUE)
            stop("Can not coerce to spam format with add.zero.rows==FALSE")
        }
        else {
            # put a hard zero in the first column of each missing row
            obj$ind <- rbind(obj$ind, cbind(ind.missing, rep(1, 
                N.missing)))
            obj$ra <- c(obj$ra, rep(0, N.missing))
        }
    }
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
    return(new("spam", entries = as.numeric(obj$ra), colindices = as.integer(obj$ind[, 
        2]), rowpointers = as.integer(ia), dimension = as.integer(obj$da)))
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
