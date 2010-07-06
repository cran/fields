# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"Krig.flplike" <- function(lambda, obj) {
    #  - log profile likelihood for lambda
    # See section 3.4 from Nychka  Spatial Processes as Smoothers paper.
    # for equation and derivation
    D2 <- obj$matrices$D[ obj$matrices$D>0]
    u2<-  obj$matrices$u[ obj$matrices$D>0]
    lD<- D2*lambda
    N2 <- length(D2)
    # MLE estimate of rho for fixed lambda
    rho.MLE<- (sum( (D2*(u2)**2)/(1+lD)))/N2
    #
    # ln determinant of    K + lambda*WI
    lnDetCov<- -sum( log(D2/(1 + lD)) )

    -1*(-N2/2 - log(2*pi)*(N2/2) - (N2/2)*log(rho.MLE) - (1/2) * lnDetCov)

  }
