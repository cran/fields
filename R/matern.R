"matern" <-
function (x = seq(0, 4 * range, len = 100), scale = 1, range = 1, 
    smoothness = 0.5) 
{
    if (!is.loaded(symbol.For("rkmat"))) {
        temp <- dyn.load(paste(FIELDS.BIN,"fields.o", sep=""), 2)
    }
    y <- x
    length <- length(x)
    theta <- c(scale, range, smoothness)
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(theta) <- "double"
    storage.mode(length) <- "integer"
    out <- .Fortran("rkmat", theta = theta, x = x, y = y, n = length)
    list(x = out$x, y = out$y, scale = scale, range = range, 
        smoothness = smoothness)
}
