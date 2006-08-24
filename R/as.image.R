"as.image" <-
function (Z, ind = NULL, grid = NULL, x = NULL, nrow = 64, ncol = 64, 
    weights = NULL, na.rm = FALSE,nx=NULL,ny=NULL) 
{

# NOTE that throughout ind is a two column integer matrix of 
# discretized locations in image matrix. 

# set some default values for arguments
#
    if (!is.null(ind)) 
        x <- ind

    if( is.null( weights)) {
      weights<- rep( 1, length(Z))}

    if (is.null(x) & is.null(grid)) {
        grid <- list(x = 1:nrow, y = 1:ncol) }
#
# use values of nx ny if passed. 
    if( !is.null(nx)) nrow<- nx
    if( !is.null(ny)) ncol<- ny


# check for missing values in Z and  na.rm==FALSE 
    if( any( is.na(Z) ) & !na.rm ) {
           stop("missing values in Z, set na.rm=TRUE")}
#
#  if there are missing overwrite Z, x, ind and weights. 
#
      if (any( is.na(Z)) ) {
         Z.good <- !is.na(Z)
         Z <- Z[Z.good]
         x <- x[Z.good, ]
         ind<- ind[Z.good,]
         weights <- weights[Z.good]
      }
#
# check for x or weights having missing. 
     if( any(is.na(weights)) | any( is.na( c( x)))) { 
           stop("missing values in weights or x")}

# if  grid is missing discretize. 
    if (!is.null(x) & is.null(grid)) {
        temp <- Krig.discretize(x, nrow, ncol)
        grid <- temp$grid
        ind <- temp$index
    }

# if both x and grid are passed discretize. 

    if (!is.null(x) & !is.null(grid)) {
        temp <- Krig.discretize(x, grid = grid)
        ind <- temp$index
    }
#
#

# m and n are the same as nrow and ncol

    m <- length(grid$x)
    n <- length(grid$y)

# find unique pixels 
#
    rep.info <- cat.matrix(ind)
    uniquerows <- !duplicated(rep.info)

#
# compute weighted means where there are replicates 
#

    if (sum(uniquerows) < length(Z)) {
        ind <- ind[uniquerows, ]
        temp <- fast.1way(rep.info, Z, w = weights)
        Z <- temp$means
        Ncell <- temp$n
        temp2 <- matrix(0, nrow = m, ncol = n)
        temp2[ind] <- Ncell
        temp3 <- matrix(NA, nrow = m, ncol = n)
        temp3[ind] <- temp$w.means
    }
    else {
        temp2 <- matrix(0, nrow = m, ncol = n)
        temp2[ind] <- 1
        temp3 <- matrix(NA, nrow = m, ncol = n)
        temp3[ind] <- 1
    }

    temp <- matrix(NA, nrow = m, ncol = n)
    temp[ind] <- Z
    call <- match.call()
    list(x = grid$x, y = grid$y, z = temp, call = call, ind = ind, 
        N = temp2, weights = temp3)
}

