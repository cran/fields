# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
 "predict.surface.derivative"<- function( object,...){
            UseMethod("predict.surface.derivative")}

 "predict.surface.derivative.default" <- function(object, grid.list = NULL, 
    nx = 80, ny = 80, ...) {
    # NOTE: this is a lean function for derivatives on 2-d grids for 2-d surfaces.
    # use predict.derivative to predict on more general problems.
  
    # without grid.list
    # default is 80X80 grid 
 
    if (is.null(grid.list) ) {
        if (is.null(object$x)) {
            stop("Need a an X matrix in the output object")
        }
        grid.list <- fields.x.to.grid(object$x, nx = nx, ny = ny, 
            xy = c(1,2))
    }
    #
    # simple check that grid list matches data
  
    if (length(grid.list)!=2) {
      stop("Only support for 2 d surfaces")}
    xg<- grid.list[[1]]
    yg<- grid.list[[2]]
    nx<- length( xg)
    ny<- length( yg)
    
    # output  'z' array to hold partials 
    out <- array(NA, c(nx,ny,2))
    #
    # explicitly loop through  grid row by row to reduce memory
    # this the position of the x grid in the list
    # fill out row by row
     for (i in 1:nx) {
#      xtemp<-cbind( rep(xg[i],ny), yg)
#      print( xtemp)
      temp <-  predict.derivative(object, x=cbind( rep(xg[i],ny), yg),  ...)
#      print( temp)
      out[i,,1] <-temp[,1]
      out[i,,2] <-temp[,2]
      
    }

   return( list(x = xg, y = yg, z = out) )

}
