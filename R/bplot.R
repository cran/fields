# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"bplot" <-
function (x, by, style = "tukey", outlier = TRUE, plot = TRUE, 
    ...) 
{
    obj <- stats.bplot(x, style = style, outlier = outlier, by = by)
    if (plot) {
        bplot.obj(obj, outlier=outlier,...)
    }
    else {
        return(obj)
    }
    invisible()
}

