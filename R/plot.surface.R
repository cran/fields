# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"plot.surface" <- function(x, main = NULL, type = "C", 
    zlab = NULL, xlab = NULL, ylab = NULL, levels = NULL, zlim = NULL, 
    graphics.reset = NULL, labcex = 0.6, add.legend = TRUE, ...) {
    obj <- x
    old.par <- par(no.readonly = TRUE)
    if (is.na(match(type, c("b", "C", "I", "p")))) {
        stop("plot type does not match b, C, I, or p.")
    }
    if (is.null(zlim)) {
        zlim = range(obj$z, na.rm = TRUE)
    }
    if (is.null(graphics.reset) & (type == "b")) {
        graphics.reset <- TRUE
    }
    else {
        graphics.reset <- FALSE
    }
    if (graphics.reset) {
        on.exit(par(old.par))
    }
    if (is.null(xlab)) {
        if (is.null(obj$xlab)) 
            xlab <- "X"
        else xlab <- obj$xlab
    }
    if (is.null(ylab)) {
        if (is.null(obj$ylab)) 
            ylab <- "Y"
        else ylab <- obj$ylab
    }
    if (is.null(zlab)) {
        if (is.null(obj$zlab)) 
            zlab <- "Z"
        else zlab <- obj$zlab
    }
    if (is.null(main)) 
        if (!is.null(obj$main)) 
            main <- obj$main
    if (type == "b") 
        set.panel(1, 2, TRUE)
    if (type == "p" | type == "b") {
        if (type == "b") {
            add.legend <- FALSE
            old.mar <- par()$mar
            par(mar = c(0, 5, 0, 0))
        }
        drape.plot(obj, xlab = xlab, ylab = ylab, zlab = zlab, 
            zlim = zlim, add.legend = add.legend, ...)
        if (!is.null(main)) 
            title(main)
    }
    if (type == "I") {
        image.plot(obj$x, obj$y, obj$z, xlab = xlab, ylab = ylab, 
            zlim = zlim, ...)
        if ((!is.null(main)) & type != "b") 
            title(main)
    }
    if (type == "b" | type == "C") {
        if (type == "b") {
            par(mar = old.mar)
        }
        image.plot(obj$x, obj$y, obj$z, xlab = xlab, ylab = ylab, 
            graphics.reset = graphics.reset, zlim = zlim, ...)
        if (is.null(levels)) 
            levels <- pretty(obj$z[!is.na(obj$z)], 5)
        contour(obj$x, obj$y, obj$z, add = TRUE, levels = levels, 
            labcex = labcex, col = "black", lwd = 2)
        if ((!is.null(main)) & type != "b") 
            title(main)
    }
    invisible()
}