# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
dyadic.check <- function(n, cut.p = 2) {
    # checks that n is of the form
    # n=p*2^m where p <= cut.p
    n2 <- as.integer(n)
    while (n2 > cut.p) {
        if (n2%%2 != 0) {
            cat(n, "must equal p*2^m where p is less than or equal to ", 
                cut.p, fill = TRUE)
            return(FALSE)
        }
        n2 <- n2/2
    }
    return(TRUE)
}
