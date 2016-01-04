# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"ceiling2" <- function(m) {
    if (m < 1) 
        return(NA)
    M <- 1
    while (M < m) {
        M <- M * 2
    }
    M
}
