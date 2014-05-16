# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"surface.Krig" <- function(object, grid.list = NULL, 
    extrap = FALSE, graphics.reset = NULL, xlab = NULL, ylab = NULL, 
    main = NULL, zlab = NULL, zlim = NULL, levels = NULL, type = "C", 
    nx = 80, ny = 80, ...) {
    ## modified so that you can give main, and ylab as arguments
    ## in ... and have them passed correctly
    out.p <- predictSurface(object, grid.list = grid.list, extrap = extrap, 
        nx = nx, ny = ny, drop.Z = TRUE)
    if (!is.null(ylab)) 
        out.p$ylab <- ylab
    if (!is.null(xlab)) 
        out.p$xlab <- xlab
    if (!is.null(zlab)) 
        out.p$zlab <- zlab
    if (!is.null(main)) 
        out.p$main <- main
    ##    else
    ##      out.p$main <- NULL
    plot.surface(out.p, type = type, graphics.reset = graphics.reset, 
        levels = levels, zlim = zlim, ...)
    invisible()
}
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"surface" <- function(object, ...) {
    UseMethod("surface")
}
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"surface.default" <- function(object, ...) {
    plot.surface(object, ...)
}
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"surface.mKrig" <- function(object, grid.list = NULL, 
    extrap = FALSE, graphics.reset = NULL, xlab = NULL, ylab = NULL, 
    main = NULL, zlab = NULL, zlim = NULL, levels = NULL, type = "C", 
    nx = 80, ny = 80, ...) {
    ## modified so that you can give main, and ylab as arguments
    ## in ... and have them passed correctly
    out.p <- predictSurface(object, grid.list = grid.list, extrap = extrap, 
        nx = nx, ny = ny, drop.Z = TRUE)
    if (!is.null(ylab)) 
        out.p$ylab <- ylab
    if (!is.null(xlab)) 
        out.p$xlab <- xlab
    if (!is.null(zlab)) 
        out.p$zlab <- zlab
    if (!is.null(main)) 
        out.p$main <- main
    ##    else
    ##      out.p$main <- NULL
    plot.surface(out.p, type = type, graphics.reset = graphics.reset, 
        levels = levels, zlim = zlim, ...)
    invisible()
}
# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"surface.surface" <- function(object, ...) {
    #
    plot.surface(object, ...)
}
