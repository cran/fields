# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"qr.q2ty" <- function(qr, y) {
    if (!is.matrix(y)) {
        y <- as.matrix(y)
    }
    dy <- dim(y)
    dq <- dim(qr$qr)
    rank <- qr$rank
    if (dy[1] != dq[1]) 
        stop("y and qr$qr should have same number of rows")
    qr.qty(qr, y)[(rank + 1):dy[1], ]
}
