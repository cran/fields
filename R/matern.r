"matern" <-
function (x = seq(0, 4 * range, , 100), scale = 1, range = 1, 
    smoothness = 0.5) 
{
    y <- x
    n <- length(x)
    theta <- c(scale, range, smoothness)
    storage.mode(x) <- "double"
    storage.mode(y) <- "double"
    storage.mode(theta) <- "double"
    storage.mode(n) <- "integer"

itest<- trunc(  smoothness)+1

ncal<- as.integer( rep( 0, n))

    out <- .Fortran("rkmat", theta = theta, x = x, y = y, n = n, 
ncal= ncal)

# set to NA any calculation that failed in the call to rkmat
out$y[ out$ncal!=itest] <- NA

 list(x = out$x, y = out$y, scale = scale, range = range,
 smoothness = smoothness)
}
