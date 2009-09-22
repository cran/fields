# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"Krig.fgcv.model" <- function(lam, obj) {
    lD <- obj$matrices$D * lam
    MSE <- sum(((obj$matrices$u * lD)/(1 + lD))^2)/length(lD)
    trA <- sum(1/(1 + lD))
    den <- (1 - (obj$cost * (trA - obj$nt - obj$offset) + obj$nt)/length(lD))
    ifelse(den > 0, obj$shat.pure.error^2 + MSE/den^2, NA)
}
