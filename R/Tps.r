"Tps" <-
function (x, Y, m = NULL, p = NULL, decomp = "WBW", scale.type = "range", ...) 
{
    x <- as.matrix(x)
    d <- ncol(x)
    if (is.null(p)) {
        if (is.null(m)) {
            m <- max(c(2, ceiling(d/2 + 0.1)))
        }
        p <- (2 * m - d)
        if (p <= 0) {
            stop(" m is too small  you must have 2*m -d >0")
        }
    }
    if (!is.null(list(...)$knots)) 
        decomp <- "DR"
    Tpscall <- match.call()
    Tpscall$cov.function <- "Thin plate spline radial basis functions (rad.cov) "
    Krig(x, Y, cov.function = rad.cov, m = m, decomp = decomp, 
        scale.type = scale.type, outputcall = Tpscall, p = p, 
        ...)
}
