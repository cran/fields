# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"predict.krig.image" <-
function (object, x, ...) 
{
    obj<- object # hack S3
    if (missing(x)) {
        temp1 <- c(obj$cov.function(obj$indexM, obj$indexM, Y = obj$delta, 
            cov.obj = obj$cov.obj))
        temp1 <- temp1[obj$rep.info]
        temp2 <- obj$beta[1] + obj$x %*% obj$beta[2:3]
        return(c(temp1 + temp2))
    }
    else {
        predict.interp.surface(obj$surface, x)
    }
}
