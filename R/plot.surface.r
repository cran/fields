"plot.surface" <-
function (obj, main = NULL, type = "b", zlab = NULL, xlab = NULL, 
    ylab = NULL, levels = NULL, zlim = NULL, graphics.reset = NULL, 
    ...) 
{
    old.par <- par(no.readonly=T)
    if (is.null(graphics.reset)&( type == "b"))
{ graphics.reset<- T}
else{ graphics.reset<- F}
    
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
        set.panel(2, 1, T)
    if (type == "p" | type == "b") {
        if (is.null(zlim)) 
            persp(obj$x, obj$y, obj$z, xlab = xlab, ylab = ylab, 
                zlab = zlab, ...)
        else persp(obj, xlab = xlab, ylab = ylab, zlab = zlab, 
            zlim = zlim, ...)
        if (!is.null(main)) 
            title(main)
    }
    if (type == "c" | type == "b") {
        if (is.null(levels)) 
            levels <- pretty(obj$z[ !is.na(obj$z)] , 5)
        contour(obj$x, obj$y, obj$z, , xlab = xlab, ylab = ylab, 
            levels = levels, ...)
    }
    if (type == "I") {
        image.plot(obj$x, obj$y, obj$z, xlab = xlab, ylab = ylab)
        if ((!is.null(main)) & type != "b") 
            title(main)
    }
    if (type == "C") {
        image.plot(obj$x, obj$y, obj$z, xlab = xlab, ylab = ylab, 
            graphics.reset = graphics.reset)
        if (is.null(levels)) 
            levels <- pretty(obj$z[ !is.na(obj$z)] , 5)
        contour(obj$x, obj$y, obj$z, add = T, levels = levels, 
            ...)
        if ((!is.null(main)) & type != "b") 
            title(main)
    }
    invisible()
}
