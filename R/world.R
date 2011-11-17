# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"world" <- function(ylim = c(-90, 90), xlim = NULL, 
    add = FALSE, asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
    eps = 0.1, col = 1, shift = FALSE, fill = FALSE, col.water = "white", 
    col.land = "darkgrey",alpha=NA, ...) {
    #
    # load world data set of land outlines
    # (should not reload if it is already loaded)
    data(world.dat)
    # check some options
    if (shift & fill) {
        stop("filling in land not implemented with shift option")
    }
    # reset colors with alpha value
   
    if (shift) {
        ind1 <- !is.na(world.dat$x)
        ind2 <- (world.dat$x < 0)
        # shift coordinates
        world.dat$x[ind2 & ind1] <- world.dat$x[ind2 & ind1] + 
            360
        # 'pick up pen' at the new edges by inserting NA's
        world.dat$x[(world.dat$x <= eps | world.dat$x >= (360 - 
            eps)) & ind1] <- NA
    }
    if (is.null(xlim)) {
        if (shift) {
            xlim <- c(0, 360)
        }
        else {
            xlim <- c(-180, 180)
        }
    }
    # create new plotting region if add is FALSE
    if (!add) {
        plot(world.dat, ylim = ylim, xlim = xlim, type = "n", 
            xaxt = xaxt, yaxt = yaxt, xlab = xlab, ylab = ylab, 
            asp = asp, ...)
    }
    # decide whether to draw lines or fill in land masses and lakes.
    if (!fill) {
        lines(world.dat, err = -1, col = col, ...)
    }
    else {
        world.color(world.dat, col.water = col.water, col.land = col.land, 
            ...)
    }
    invisible()
}
world.color <- function(obj,xlim = c(-180, 180), ylim = c(-90, 
    90), col.water = "white", col.land = "darkgrey", ...) {

    # obj should be world.dat! This is done to avoid a funny 
    # warning in the R package checking.

    # logicals for land and lakes in case these need to
    # be modified.
    land <- TRUE
    lakes <- TRUE
    # first add ocean color as background
    # find the size of current region but don't let it get bigger than
    #lon lat range
    #if(!no.ocean){
    #rect(xlim[1], ylim[1], xlim[2], ylim[2], col = col.water, 
    #    border = NA)}
    
    # find separate polygons of land masses and islands in world.dat
    # these are indicated by NA's
    ind <- (1:length(obj$x))[is.na(obj$x)]
    ind <- c(1, ind)
    N <- length(ind) - 1
    lakes.id <- c(46, 53, 25, 26, 28, 27, 4, 47, 48, 51, 49, 
        50)
    land.id <- (1:N)[-lakes.id]
    # loop through polygons of land
    if (land) {
        for (k in land.id) {
            tempi <- ind[k]:ind[k + 1]
            polygon(list(x = obj$x[tempi], y = obj$y[tempi]), 
                col = col.land, border = col.land, ...)
        }
    }
    #
    # add in interior of antarctica as a
    # rectangle below farthest south coastline
    ytemp <- obj$y[ind[21] + 1]
    rect(-180, -90, 180, ytemp, col = col.land, border = col.land, 
        ...)
    # loop through polygons of lakes do this second so
    # color over writes the land
    if (lakes) {
        for (k in lakes.id) {
            tempi <- ind[k]:ind[k + 1]
            polygon(list(x = obj$x[tempi], y = obj$y[tempi]), 
                col = col.water, border = col.water, ...)
        }
    }
}
world.land <- function( col.water = "white", col.land = "darkgrey",alpha=NA, ...) {
    data(world.dat)
    obj<- world.dat
     if( !is.na( alpha) ){
    #    
      col.water<- rgb( t(col2rgb(col.water, TRUE)), alpha=255*alpha, maxColorValue=255)
      col.land<-  rgb( t(col2rgb(col.land, TRUE)), alpha=255*alpha, maxColorValue=255)
    }
    # find separate polygons of land masses and islands in world.dat
    # these are indicated by NA's
    ind <- (1:length(obj$x))[is.na(obj$x)]
    ind <- c(1, ind)
    N <- length(ind) - 1
    lakes.id <- c(46, 53, 25, 26, 28, 27, 4, 47, 48, 51, 49, 
        50)
    land.id <- (1:N)[-lakes.id]
    # loop through polygons of land
   
        for (k in land.id) {
            tempi <- ind[k]:ind[k + 1]
            polygon(list(x = obj$x[tempi], y = obj$y[tempi]), 
                col = col.land, border = col.land, ...)
        }
    
    #
    # add in interior of antarctica as a
    # rectangle below farthest south coastline
    ytemp <- obj$y[ind[21] + 1]
    rect(-180, -90, 180, ytemp, col = col.land, border = col.land, 
        ...)
    # loop through polygons of lakes do this second so
    # color over writes the land
   
        for (k in lakes.id) {
            tempi <- ind[k]:ind[k + 1]
            polygon(list(x = obj$x[tempi], y = obj$y[tempi]), 
                col = col.water, border = col.water, ...)
        }
    
}

in.land.grid<- function (grid.list)
{
data(world.dat)
ind0<-is.na(world.dat$x)|is.na(world.dat$x)
 ind <- (1:length(world.dat$x))[ind0]
    ind <- c(0, ind)
    N <- length(ind) - 1
    lakes.id <- c(46, 53, 25, 26, 28, 27, 4, 47, 48, 51, 49, 
        50)
    start.id<- ind[1:N] +1
    end.id<-   ind[2:(N+1)] -1

    land.id <- (1:N)[-lakes.id]
        m<- length( grid.list$x)
        n<- length( grid.list$y)
        land.mask<- matrix( FALSE, m,n)
        for (k in land.id) {
            tempi <- start.id[k] : end.id[k]
            xp<- cbind(world.dat$x[tempi],  world.dat$y[tempi] )
            land.mask<- land.mask | in.poly.grid( grid.list,xp, convex.hull=FALSE)  
         }

        for (k in lakes.id) {
            tempi <- start.id[k] : end.id[k]
            xp<- cbind(world.dat$x[tempi],  world.dat$y[tempi] )
            land.mask<- land.mask & !(in.poly.grid( grid.list,xp))
        }
 return( land.mask)
}


