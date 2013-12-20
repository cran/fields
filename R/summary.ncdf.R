# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"summary.ncdf" <- function (object, ...) 
{
     tempList<- NULL
     varNames<- NULL
#
    cat("DIMENSIONS", fill=TRUE)
    for (i in names(object$dim) ) {
        vname = i
        ndims = length(object$dim[[i]]$vals)
         cat(vname, " has size", ndims, fill=TRUE)
    }
    cat(fill=TRUE) 
    cat("VARIABLES", fill=TRUE)
    for (i in 1:object$nvars) {
        vname = object$var[[i]]$name
        ndims = object$var[[i]]$ndims
        dimstring = paste(vname, "( variable ",i , 
            ") has shape")
        dimTemp<- NULL
        for (j in 1:ndims) {
            dimTemp<- c( dimTemp, object$var[[i]]$dim[[j]]$len)
        }
        temp<- ( dimTemp)
        varNames<- c(varNames, vname)
        tempList<- c( tempList, list(dimTemp))
        if( is.null(dimTemp) ){
           dimTemp<- NA}
        cat( i,":",  vname, "has size ", dimTemp, sep=" ", fill = TRUE)
    }
     names(tempList) <- varNames
     invisible( tempList)
}
