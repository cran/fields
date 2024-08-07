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
"quilt.plot" <- function(x, y, z, nx = 64, ny = 64, 
     grid = NULL, add.legend = TRUE, add = FALSE, nlevel=64, 
    col = tim.colors(nlevel), nrow = NULL, ncol = NULL, FUN=NULL,
    plot=TRUE, na.rm=FALSE, boundary.grid = FALSE,  ...) {
    #
    # note that nrow and ncol refer to the resulting 'image format' for plotting.
    # here the x values are the rows and the y values are the columns
    # FUN = NULL means the weighted means are found for each grid cell
    if( !is.null(nrow)|!is.null(nrow)){
      nx<- nrow
      ny<- ncol
     }
    x <- as.matrix(x)
    if (ncol(x) == 2) {
        z <- y
    }
    if (ncol(x) == 1) {
        x <- cbind(x, y)
    }
    if (ncol(x) == 3) {
        z <- x[, 3]
        x <- x[, 1:2]
    }
    # at this point x should be a 2 column matrix of x-y locations
    #  z is a vector or one column matrix of the z values.
    #discretize data
    out.p <- as.image(z, x = x, nx = nx, ny = ny, 
        grid = grid, FUN=FUN, na.rm=na.rm, boundary.grid = boundary.grid)
    # besides the image information this list has the indices that 
    # map each z value to a grid box
    #    
    # plot it
    if( plot){
    if (add.legend) {
        image.plot(out.p, nlevel = nlevel, col = col, add = add, ...)
    }
    else {
        image(out.p, col = col, add = add, ...)
    }
    }
    invisible(out.p)
}
