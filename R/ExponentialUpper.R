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
ExponentialUpper = function(distMat, range = 1, alpha = 1/range,
                            theta=NULL) {
  # theta argument has been depreciated.
  if( !is.null( theta)){
    aRange<- theta
  }
  # Evaluates the exponential covariance function over the upper triangle of the distance matrix
  
  if(nrow(distMat) != ncol(distMat))
    stop('distance matrix is non-symmetric.  Should not be calling ExponentialUpper.')
  
  return(.Call("ExponentialUpperC", as.double(distMat), as.integer(nrow(distMat)), as.double(alpha), PACKAGE = "fields"))
  
  #convert ans to standard matrix
  #ans = ans[[1]]
  #dim(ans) = dim(distMat)
  
  #return(ans)
}
