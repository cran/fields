"bplot" <-
function (x, style = "tukey", outlier = T, plot = T, ...) 
{
    obj <- stats.bplot(x, style = style, outlier = outlier)
    if (plot) {
        bplot.obj(obj, ...)
    }
    else {
        return(obj)
    }
    invisible()
}
