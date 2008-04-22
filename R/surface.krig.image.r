# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"surface.krig.image" <-
function (obj, grid.list = NA, extrap = TRUE, graphics.reset = FALSE, 
    xlab = NULL, ylab = NULL, main = NULL, zlab = NULL, zlim = NULL, 
    levels = NULL, ptype = "I", ...) 
{
    old.par <- par("mfrow", "oma")
    if (graphics.reset) 
        on.exit(par(old.par))
    if (is.na(grid.list)) {
        out.p <- obj$surface
    }
    else {
        out.p <- predict.surface(obj, grid.list = grid.list, 
            extrap = extrap)
    }
    if (!is.null(ylab)) 
        out.p$ylab <- ylab
    if (!is.null(xlab)) 
        out.p$xlab <- xlab
    if (!is.null(zlab)) 
        out.p$zlab <- zlab
    if (!is.null(main)) 
        out.p$main <- main
    plot.surface(out.p, type = ptype, graphics.reset = graphics.reset, 
        levels = levels, zlim = zlim, ...)
    invisible()
}
