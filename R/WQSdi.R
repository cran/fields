# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WQSdi" <- function(x, transpose = FALSE) {
    if (!transpose) 
        t(WQSi(t(WQSi(x))))
    else {
        WQSi.T(t(WQSi.T(t(x))))
    }
}
