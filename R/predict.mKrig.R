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
predict.mKrig <- function(object, xnew = NULL, ynew = NULL, grid.list=NULL,
                          derivative = 0, Z = NULL, drop.Z = FALSE, just.fixed = FALSE,
                          collapseFixedEffect = object$collapseFixedEffect, 
                          ...) {
  # the main reason to pass new args to the covariance is to increase
  # the temp space size for sparse multiplications
  # other optional arguments that typically describe the covariance function 
  # from mKrig are passed along in the list object$args
  cov.args <- list(...)
  # predict at observation locations by default
  if( !is.null(grid.list)){
    xnew<- make.surface.grid(grid.list)
  }
  if (is.null(xnew)) {
    xnew <- object$x
  }
  if (is.null(Z) & (length(object$ind.drift) >0 )) {
    Z <- object$Tmatrix[, !object$ind.drift]
  }
  if (!is.null(ynew)) {
    coef.hold <- mKrig.coef(object, ynew,
                            collapseFixedEffect=collapseFixedEffect)
    c.coef <- coef.hold$c.coef
    beta <- coef.hold$beta
  }
  else {
    c.coef <- object$c.coef
    beta <- object$beta
  }
  # fixed part of the model this a polynomial of degree m-1
  # Tmatrix <- fields.mkpoly(xnew, m=object$m)
  # only do this if nt>0, i.e. there is a fixed part.
  #
  
  if( object$nt>0){
    if (derivative == 0) {
      if (drop.Z | object$nZ == 0) {
        # just evaluate polynomial and not the Z covariate
        temp1 <- fields.mkpoly(xnew, m = object$m) %*% 
          beta[object$ind.drift, ]
      }
      else {
        if( nrow( xnew) != nrow(as.matrix(Z)) ){
          stop(paste("number of rows of covariate Z",
                     nrow(as.matrix(Z)), 
                     " is not the same as the number of locations",
                     nrow( xnew) )
          )
        }
        temp0 <-  cbind(fields.mkpoly(xnew, m = object$m),as.matrix(Z)) 
        temp1 <- temp0 %*% beta
      }
    }
    else {
      if (!drop.Z & object$nZ > 0) {
        stop("derivative not supported with Z covariate included")
      }
      temp1 <- fields.derivative.poly(xnew, m = object$m, beta[object$ind.drift, 
      ])
    }
    if (just.fixed) {
      return(temp1)
    }
  }  
  # add nonparametric part. Covariance basis functions
  # times coefficients.
  # syntax is the name of the function and then a list with
  # all the arguments. This allows for different covariance functions
  # that have been passed as their name.
  if (derivative == 0) {
    # argument list are the parameters and other options from mKrig
    #  locations anbetaficients,
    temp2 <- do.call(object$cov.function.name, c(object$args, 
                                                 list(x1 = xnew, x2 = object$knots, C = c.coef), cov.args))
  }
  else {
    temp2 <- do.call(object$cov.function.name, c(object$args, 
                                                 list(x1 = xnew, x2 = object$knots, C = c.coef, derivative = derivative), 
                                                 cov.args))
  }
  # add two parts together and coerce to vector
  if( object$nt>0){
    return((temp1 + temp2))
  }
  else{
    return(  temp2)
  }
}
