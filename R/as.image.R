# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"as.image" <- function(Z, ind = NULL, grid = NULL, 
    x = NULL,  weights = rep(1, length(Z)), na.rm = FALSE, 
    nx = 64, ny = 64, boundary.grid = FALSE, nrow = NULL, ncol = NULL) {
    # NOTE that throughout ind is a two column integer matrix of
    # discretized locations in the image matrix.
    # Thanks to J. Rougier for fixing bugs in this function.
    # set some default values for arguments
    #
    # coerce Z to a vector
    Z <- c(Z)
    if( !is.null(ind)){
      x<- ind
     }
    if( !is.null(nrow)&!is.null(ncol)){
      nx<- nrow
      ny<- ncol
    }
    #
    # check for x or weights having missing values
    # we do not like these ...
    if (any(is.na(weights)) | any(is.na(c(x)))) {
        stop("missing values in weights or x")
    }
    # discretize locations to grid boxes
    # this function will also create a default grid based on range of
    # locations if grid is NULL
    #
    temp <- discretize.image(x, m = nx, n = ny, grid = grid, 
        boundary.grid = boundary.grid)
    grid <- temp$grid
    # index is a two component list that indexes  the x and y grid points.
    # points outside of grid are assigned as NA
    #
    # empty image matrices to hold weights and  weighted means
     w<- z <- matrix( NA, nrow=temp$m, ncol=temp$n)
     # find stats
     tempz<- tapply( Z*weights, temp$index, sum, na.rm=FALSE )
     tempw<- tapply( weights, temp$index, sum, na.rm=FALSE)
     # these are the indices that are represented by the locations
     # they may not include the entire set ( 1:nx and 1:ny)
     # so define what they do have.
  
     # insert the tabled values into the right rows and columns.
      z[ temp$ix, temp$iy] <- tempz/ tempw
      w[ temp$ix, temp$iy] <- tempw
     # save call
     # xd created because it is a  pain to do otherwise and handy to have
    call <- match.call()
    list(x = grid$x, y = grid$y, z = z, call = call, ind = cbind(temp$index[[1]], temp$index[[2]]) , 
        weights = w, xd = cbind(grid$x[temp$index[[1]]], grid$y[temp$index[[2]]] ), 
        call = match.call() )
}
