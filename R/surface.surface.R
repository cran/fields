"surface.surface" <-
function (obj, lab = NA, type = "b", zlab, xlab, ylab, graphics.reset = FALSE, 
    ...) 
{
    old.par <- par("mfrow", "oma")
    if (graphics.reset) {
        on.exit(par(old.par))
        par(xpd = TRUE)
    }
    if (is.null(obj$xlab)) 
        obj$xlab <- "X"
    if (is.null(obj$ylab)) 
        obj$ylab <- "Y"
    if (missing(zlab)) {
        zlab <- "Z"
    }
    if (!missing(xlab)) {
        obj$xlab <- xlab
    }
    if (!missing(ylab)) {
        obj$xlab <- ylab
    }
    if (is.na(lab) & !is.null(obj$main)) {
        lab <- obj$main
    }
    if (type == "b") 
        set.panel(2, 1, TRUE)
    if (type == "p" | type == "b") {
        persp(obj, xlab = obj$xlab, ylab = obj$ylab, zlab = zlab, 
            ...)
        if (!is.na(lab)) {
            title(lab)
        }
    }
    if (type == "c" | type == "b") {
        contour(obj, xlab = obj$xlab, ylab = obj$ylab, ...)
        if (!is.na(lab) & type != "b") {
            title(lab)
        }
    }
    if (type == "i" | type == "I") {
        image.plot(x = obj$x, y = obj$y, z = obj$z, xlab = obj$xlab, 
            ylab = obj$ylab, ...)
        if (!is.na(lab) & type != "b") {
            title(lab)
        }
        if (type == "I") {
            contour(obj, add = TRUE)
        }
    }
    invisible()
}
