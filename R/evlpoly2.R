# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

`fields.evlpoly2` <-
function (x, coef, ptab) 
{

# evaluates polynomial at x values with coefficients coef[i] and powers i-1
# 
    n <- nrow(x)
    nd<- ncol( x)
     
    J<- nrow( ptab)
     if( length( coef) != J) { 
        stop("coefficients not same length as ptab rows")}    

    results<- rep(0, n)

    temp<- .Fortran("evlpoly2", 
            x=as.double(x), n=as.integer(n), nd=as.integer( nd),
            ptab=as.integer( ptab),  j=as.integer(J),
            coef=as.double(coef),
            results=as.double(results)
           , PACKAGE = "fields")$results
    return( temp)
}

