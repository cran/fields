"predict.qsreg" <-
function (out, x, derivative = 0, model = out$ind.cv.ps) 
{
    if (missing(x)) 
        x <- out$x
    c(splint(out$predicted$x, out$predicted$y[, model], x, derivative = derivative))
}
