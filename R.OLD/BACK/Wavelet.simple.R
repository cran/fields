# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
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
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WD42d" <- function(x) {
    t(WD4(t(WD4(x))))
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WD42di" <- function(x) {
    WD4i(t(WD4i(t(x))))
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WD4i" <- function(x) {
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    n <- nrow(x)
    m <- ncol(x)
    X <- matrix(0, nrow = n, ncol = m)
    X[seq(1, n, 2), ] <- x[1:(n/2), ]
    X[seq(2, n, 2), ] <- x[(n/2 + 1):n, ]
    x <- X
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
    D4tr.L <- matrix(0, nrow = 3, ncol = 4)
    D4tr.L[1, ] <- c(D4.le[1, 1], D4.le[2, 1], D4.le[3, 1], D4.le[4, 
        1])
    D4tr.L[2, ] <- c(D4.le[1, 2], D4.le[2, 2], D4.le[3, 2], D4.le[4, 
        2])
    D4tr.L[3, ] <- c(D4.le[1, 3], D4.le[2, 3], D4.le[3, 3], D4.le[4, 
        3])
    D4tr.R <- matrix(0, nrow = 3, ncol = 4)
    D4tr.R[1, ] <- c(D4.re[1, 3], D4.re[2, 3], D4.re[3, 3], D4.re[4, 
        3])
    D4tr.R[2, ] <- c(D4.re[1, 4], D4.re[2, 4], D4.re[3, 4], D4.re[4, 
        4])
    D4tr.R[3, ] <- c(D4.re[1, 5], D4.re[2, 5], D4.re[3, 5], D4.re[4, 
        5])
    if (n == 8) {
        D4tr.t <- matrix(0, nrow = 2, ncol = 4)
        D4tr.t[1, ] <- c(D4.le[3, 4], D4.le[4, 4], D4.re[1, 1], 
            D4.re[2, 1])
        D4tr.t[2, ] <- c(D4.le[3, 5], D4.le[4, 5], D4.re[1, 2], 
            D4.re[2, 2])
    }
    else {
        D4tr.t <- matrix(0, 4, 4)
        D4tr.t[1, ] <- c(D4.le[3, 4], D4.le[4, 4], D4[1], D4.r[1])
        D4tr.t[2, ] <- c(D4.le[3, 5], D4.le[4, 5], D4[2], D4.r[2])
        D4tr.t[3, ] <- c(D4[3], D4.r[3], D4.re[1, 1], D4.re[2, 
            1])
        D4tr.t[4, ] <- c(D4[4], D4.r[4], D4.re[1, 2], D4.re[2, 
            2])
    }
    iD4 <- c(D4[3], D4.r[3], D4[1], D4.r[1])
    iD4.r <- c(D4[4], D4.r[4], D4[2], D4.r[2])
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[1:3, ] <- D4tr.L %*% x[1:4, ]
    tmp[(n - 2):n, ] <- D4tr.R %*% x[(n - 3):n, ]
    if (n == 8) {
        tmp[4:5, ] <- D4tr.t %*% x[3:6, ]
    }
    else {
        tmp[4:5, ] <- D4tr.t[1:2, ] %*% x[3:6, ]
        tmp[(n - 4):(n - 3), ] <- D4tr.t[3:4, ] %*% x[(n - 5):(n - 
            2), ]
        indx <- seq(6, (n - 6), 2)
        tmp[indx, ] <- iD4[1] * x[indx - 1, ] + iD4[2] * x[indx, 
            ] + iD4[3] * x[indx + 1, ] + iD4[4] * x[indx + 2, 
            ]
        tmp[indx + 1, ] <- iD4.r[1] * x[indx - 1, ] + iD4.r[2] * 
            x[indx, ] + iD4.r[3] * x[indx + 1, ] + iD4.r[4] * 
            x[indx + 2, ]
    }
    iwt <- tmp
    iwt
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WQS" <- function(x) {
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D.smooth <- c(-1, 3, 3, -1)
    D.le <- matrix(c(2, 3, -1, 2, -3, 1), nrow = 2, ncol = 3, 
        byrow = TRUE)
    D.re <- matrix(c(-1, 3, 2, -1, 3, -2), nrow = 2, ncol = 3, 
        byrow = TRUE)
    D.rough <- c(-1, 3, -3, 1)
    n <- dim(x)[1]
    m <- dim(x)[2]
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[1:2, ] <- D.le %*% x[1:3, ]
    tmp[(n - 1):n, ] <- D.re %*% x[(n - 2):n, ]
    stuff <- n - 4
    indx <- seq(2, stuff, 2)
    tmp[indx + 1, ] <- D.smooth[1] * x[indx, ] + D.smooth[2] * 
        x[indx + 1, ] + D.smooth[3] * x[indx + 2, ] + D.smooth[4] * 
        x[indx + 3, ]
    tmp[indx + 2, ] <- D.rough[1] * x[indx, ] + D.rough[2] * 
        x[indx + 1, ] + D.rough[3] * x[indx + 2, ] + D.rough[4] * 
        x[indx + 3, ]
    tmp[c(seq(1, n, 2), seq(2, n, 2)), ]
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WQS.T" <- function(x) {
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D.smooth <- c(3, -3, -1, -1)
    D.le <- matrix(c(2, 2), nrow = 1, ncol = 2, byrow = TRUE)
    D.re <- matrix(c(2, -2), nrow = 1, ncol = 2, byrow = TRUE)
    D.rough <- c(-1, 1, 3, 3)
    n <- dim(x)[1]
    m <- dim(x)[2]
    x[c(seq(1, n, 2), seq(2, n, 2)), ] <- x
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[1, ] <- D.le %*% x[1:2, ]
    tmp[n, ] <- D.re %*% x[(n - 1):n, ]
    stuff <- n - 3
    indx <- seq(1, stuff, 2)
    tmp[indx + 1, ] <- D.smooth[1] * x[indx, ] + D.smooth[2] * 
        x[indx + 1, ] + D.smooth[3] * x[indx + 2, ] + D.smooth[4] * 
        x[indx + 3, ]
    tmp[indx + 2, ] <- D.rough[1] * x[indx, ] + D.rough[2] * 
        x[indx + 1, ] + D.rough[3] * x[indx + 2, ] + D.rough[4] * 
        x[indx + 3, ]
    tmp
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WQS.basis" <- function(N, cut.n = 8) {
    Wtransform(diag(1, N), inv = TRUE, cut.min = cut.n)
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WQS2d" <- function(x, transpose = FALSE) {
    if (!transpose) 
        t(WQS(t(WQS(x))))
    else {
        WQS.T(t(WQS.T(t(x))))
    }
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WQS2di" <- function(x, transpose = FALSE) {
    if (!transpose) 
        t(WQSi(t(WQSi(x))))
    else {
        WQSi.T(t(WQSi.T(t(x))))
    }
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WQSdi" <- function(x, transpose = FALSE) {
    if (!transpose) 
        t(WQSi(t(WQSi(x))))
    else {
        WQSi.T(t(WQSi.T(t(x))))
    }
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WQSi" <- function(x) {
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D.smooth <- c(3, -3, 1, 1)/16
    D.le <- matrix(c(4, 4), nrow = 1, ncol = 2, byrow = TRUE)/16
    D.re <- matrix(c(4, -4), nrow = 1, ncol = 2, byrow = TRUE)/16
    D.rough <- c(1, -1, 3, 3)/16
    n <- dim(x)[1]
    m <- dim(x)[2]
    x[c(seq(1, n, 2), seq(2, n, 2)), ] <- x
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[1, ] <- D.le %*% x[1:2, ]
    tmp[n, ] <- D.re %*% x[(n - 1):n, ]
    stuff <- n - 3
    indx <- seq(1, stuff, 2)
    tmp[indx + 1, ] <- D.smooth[1] * x[indx, ] + D.smooth[2] * 
        x[indx + 1, ] + D.smooth[3] * x[indx + 2, ] + D.smooth[4] * 
        x[indx + 3, ]
    tmp[indx + 2, ] <- D.rough[1] * x[indx, ] + D.rough[2] * 
        x[indx + 1, ] + D.rough[3] * x[indx + 2, ] + D.rough[4] * 
        x[indx + 3, ]
    tmp
}
# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"WQSi.T" <- function(x) {
    if (!is.matrix(x)) 
        x <- matrix(x, nrow = length(x), ncol = 1)
    D.smooth <- c(1, 3, 3, 1)/16
    D.le <- matrix(c(4, 3, 1, 4, -3, -1), nrow = 2, ncol = 3, 
        byrow = TRUE)/16
    D.re <- matrix(c(1, 3, 4, 1, 3, -4), nrow = 2, ncol = 3, 
        byrow = TRUE)/16
    D.rough <- c(1, 3, -3, -1)/16
    n <- dim(x)[1]
    m <- dim(x)[2]
    tmp <- matrix(NA, nrow = n, ncol = m)
    tmp[1:2, ] <- D.le %*% x[1:3, ]
    tmp[(n - 1):n, ] <- D.re %*% x[(n - 2):n, ]
    stuff <- n - 4
    indx <- seq(2, stuff, 2)
    tmp[indx + 1, ] <- D.smooth[1] * x[indx, ] + D.smooth[2] * 
        x[indx + 1, ] + D.smooth[3] * x[indx + 2, ] + D.smooth[4] * 
        x[indx + 3, ]
    tmp[indx + 2, ] <- D.rough[1] * x[indx, ] + D.rough[2] * 
        x[indx + 1, ] + D.rough[3] * x[indx + 2, ] + D.rough[4] * 
        x[indx + 3, ]
    tmp[c(seq(1, n, 2), seq(2, n, 2)), ]
}
