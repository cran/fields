#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2024 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.com,
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2
##END HEADER
"stationary.taper.cov" <- function(x1, x2=NULL, Covariance = "Exponential", 
    Taper = "Wendland", Dist.args = NULL, Taper.args = NULL, 
    aRange = 1, V = NULL, C = NA, marginal = FALSE, spam.format = TRUE, 
    verbose = FALSE, theta=NULL, ...) {
     # theta argument has been deopreciated.
  if( !is.null( theta)){
    aRange<- theta
  }
                                        # get covariance function arguments from call
    
    Cov.args <- list(...)
    # coerce x1 and x2 to matrices
    if (is.data.frame(x1)) 
        x1 <- as.matrix(x1)
    if (!is.matrix(x1)) 
        x1 <- matrix(c(x1), ncol = 1)
    if (is.null(x2)) 
        x2 <- x1
    if (is.data.frame(x2)) 
        x2 <- as.matrix(x1)
    if (!is.matrix(x2)) 
        x2 <- matrix(c(x2), ncol = 1)
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    # Default taper arguments that are particular to the Wendland.
    # Make sure dimension argument is added.
    if (Taper == "Wendland") {
        if (is.null(Taper.args)) {
            Taper.args <- list(aRange = 1, k = 2, dimension = ncol(x1))
        }
        if (is.null(Taper.args$dimension)) {
            Taper.args$dimension <- ncol(x1)
        }
    }
    #
    # Add in general defaults for taper arguments if not Wendland
    #  aRange = 1.0 is the default range for the taper.
    if (is.null(Taper.args)) {
        Taper.args <- list(aRange = 1)
    }
    #
    # separate out a single scalar transformation and a
    # more complicated scaling and rotation.
    # this is done partly to have the use of great circle distance make sense
    # by applying the scaling  _after_ finding the distance.
    #
    # flag for great circle distance
    great.circle <- ifelse(is.null(Dist.args$method), FALSE, 
        Dist.args$method == "greatcircle")
    # check form of aRange
    if (length(aRange) > 1) {
        stop("aRange as a matrix has been depreciated,  use the V argument")
    }
    #
    # following now treats V as a full matrix for scaling and rotation.
    #
    if (!is.null(V)) {
        # try to catch error of mixing great circle distance with a
        # linear scaling of coordinates.
        if (aRange != 1) {
            stop("can't specify both aRange and V!")
        }
        if (great.circle) {
            stop("Can not mix great circle distance\nwith general scaling (V argument or vecotr of aRange's)")
        }
        x1 <- x1 %*% t(solve(V))
        x2 <- x2 %*% t(solve(V))
    }
    #
    # locations are now scaled and rotated correctly
    # copy taper range
    if (great.circle) {
        # set the delta cutoff to be in scale of angular latitude.
        # figure out if scale is in miles or kilometers
        miles <- ifelse(is.null(Dist.args$miles), TRUE, Dist.args$miles)
        delta <- (180/pi) * Taper.args$aRange/ifelse(miles, 3963.34, 
            6378.388)
    }
    else {
        delta <- Taper.args$aRange
    }
    if (length(delta) > 1) {
        stop("taper range must be a scalar")
    }
    #NOTE tapering is applied to the _scaled_ locations.
    # now apply covariance function to pairwise distance matrix, or multiply
    # by C vector or just find marginal variance
    if (!marginal) {
        # find nearest neighbor distances based on taper threshhold.
        # This is hardwired to 'nearest.dist' function from spam.
        # note that delta is taken from the taper range not aRange or V
        sM <- do.call("nearest.dist", c(list(x1, x2, delta = delta, 
            upper = NULL), Dist.args))
        # sM@entries are the pairwise distances up to distance taper.range.
        # apply covariance and taper to these.
        # note rescaling by aRange and taper ranges.
        sM@entries <- do.call(Covariance, c(list(d = sM@entries/aRange), 
            Cov.args)) * do.call(Taper, c(list(d = sM@entries), 
            Taper.args))
        # if verbose print out each component separately
        if (verbose) {
            print(sM@entries/aRange)
            print(do.call(Covariance, c(list(d = sM@entries/aRange), 
                Cov.args)))
            print(do.call(Taper, c(list(d = sM@entries), Taper.args)))
        }
        if (is.na(C[1])) {
            # decide whether to return sM in spam sparse form or as a full matrix
            if (spam.format) {
                return(sM)
            }
            else {
                return(as.matrix(sM))
            }
        }
        else {
            # other option is to sparse multiply cross covariance by C
            return(sM %*% C)
        }
    }
    else {
        # find marginal variance and return  a vector.
        tau2 <- do.call(Covariance, c(list(d = 0), Cov.args)) * 
            do.call(Taper, c(list(d = 0), Taper.args))
        return(rep(tau2, nrow(x1)))
    }
    # should not get here!
}
