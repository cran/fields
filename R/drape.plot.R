"drape.plot" <-
function (x, y, z, z2 = NULL, col = tim.colors(64), zlim = NULL, 
    add.legend = TRUE, horizontal = TRUE, theta = 30, phi = 20, 
    ...) 
{
    if (is.list(x)) {
        z <- x$z
        y <- x$y
        x <- x$x
    }
    NC <- length(col)
    M <- nrow(z)
    N <- ncol(z)
    if (!is.null(z2)) {
        M2 <- nrow(z2)
        N2 <- ncol(z2)
        if ((M != M2) | (N != N2)) {
            stop("draping matrix dimensions must match z")
        }
    }
    else {
        z2 <- z
    }
    if (is.null(zlim)) {
        zlim <- range(c(z2), na.rm = TRUE)
    }
    zcol <- drape.color(z2, col = col, zlim = zlim)
    pm <- persp(x, y, z, theta = theta, phi = phi, col = zcol, 
        ...)
    if (add.legend) {
        image.plot(zlim = zlim, legend.only = TRUE, col = col, 
            horizontal = horizontal)
    }
    invisible(pm)
}

