`fields.evlpoly` <-
function (x, coef) 
{

# evaluates polynomial at x values with coefficients coef[i] and powers i-1
# 
    n <- length(x)
    J<- length( coef)
    results<- rep(0, n)

    temp<- .Fortran("evlpoly", x=as.double(x), n=as.integer(n), 
            coef=as.double(coef), j=as.integer(J), 
            results=as.double(results),
            PACKAGE = "fields")$results

    return( temp)
}

