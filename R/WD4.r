# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WD4" <- function(x) {
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D4 <- c(0.482962913145, 0.836516303738, 0.224143868042, -0.129409522551)
    D4.le <- matrix(c(0.6033325119, -0.7965435169, 0.03751746045, 
        0.01003722456, 0.6908955318, 0.546392714, 0.4573276599, 
        0.1223510431, -0.3983129977, -0.2587922483, 0.8500881025, 
        0.2274281117, 0, 0, 0.223820357, -0.8366029212, 0, 0, 
        -0.1292227434, 0.4830129218), nrow = 4, ncol = 5, byrow = FALSE)
    D4.r <- c(-0.129409522551, -0.224143868042, 0.836516303738, 
        -0.482962913145)
    D4.re <- matrix(c(0.4431490496, 0.231557595, 0, 0, 0.7675566693, 
        0.4010695194, 0, 0, 0.3749553316, -0.7175799994, 0.2303890438, 
        -0.5398225007, 0.1901514184, -0.3639069596, 0.434896998, 
        0.801422962, -0.1942334074, 0.3717189665, 0.8705087534, 
        -0.2575129195), nrow = 4, ncol = 5, byrow = FALSE)
    n <- dim(x)[1]
    m <- dim(x)[2]
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[1:4, ] <- D4.le %*% x[1:5, ]
    tmp[(n - 3):n, ] <- D4.re %*% x[(n - 4):n, ]
    stuff <- n - 6
    indx <- seq(4, stuff, 2)
    tmp[indx + 1, ] <- D4[1] * x[indx, ] + D4[2] * x[indx + 1, 
        ] + D4[3] * x[indx + 2, ] + D4[4] * x[indx + 3, ]
    tmp[indx + 2, ] <- D4.r[1] * x[indx, ] + D4.r[2] * x[indx + 
        1, ] + D4.r[3] * x[indx + 2, ] + D4.r[4] * x[indx + 3, 
        ]
    wt <- rbind(tmp[seq(1, n, 2), ], tmp[seq(2, n, 2), ])
    wt
}
