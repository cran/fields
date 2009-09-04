# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"minimax.crit" <- function(obj, des = TRUE, R) {
    R <- as.matrix(R)
    id <- 1:nrow(R)
    if (des) 
        Dset <- attr(obj, "best.id")
    else Dset <- obj
    Cset <- id[-Dset]
    dist.mat <- rdist(R[Cset, ], R[Dset, ])
    mM.crit <- max(apply(dist.mat, 1, min))
    mM.crit
}
