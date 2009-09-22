# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"crop.image" <- function(obj, loc = NULL, ...) {
    if (is.null(loc)) {
        image.plot(obj, ...)
        loc <- get.rectangle()
    }
    # coerce to midpoints
    m <- nrow(obj$z)
    n <- ncol(obj$z)
    nx <- length(obj$x)
    ny <- length(obj$y)
    if (nx != m) {
        obj$x <- (obj$x[1:m] + obj$x[2:(m + 1)])/2
    }
    if (ny != n) {
        obj$y <- (obj$y[1:n] + obj$x[2:(n + 1)])/2
    }
    # coerce loc to x,y list format if matrix  or data frame
    if (is.matrix(loc) | is.data.frame(loc)) {
        if (ncol(loc) != 2) {
            stop("loc must have two columns\n(for x and y coordinates )")
        }
        loc <- list(x = loc[, 1], y = loc[, 2])
    }
    x <- obj$x
    y <- obj$y
    N <- length(x)
    xr <- range(loc$x)
    xtest <- range(x)
    if (xr[1] < xtest[1] | xr[2] > xtest[2]) {
        stop("cropping outside ranges of x values")
    }
    x1 <- max((1:N)[xr[1] >= x])
    x2 <- min((1:N)[xr[2] <= x])
    N <- length(y)
    yr <- range(loc$y)
    ytest <- range(y)
    if (yr[1] < ytest[1] | yr[2] > ytest[2]) {
        stop("cropping outside ranges of y values")
    }
    y1 <- max((1:N)[yr[1] >= y])
    y2 <- min((1:N)[yr[2] <= y])
    list(x = obj$x[x1:x2], y = obj$y[y1:y2], z = obj$z[x1:x2, 
        y1:y2])
}
average.image <- function(obj, Q = 2) {
    # fast method to sum over a QXQ block in image.
    # Q is the number of elements to average over in each dimension
    # e.g.  Q=5 --  blocks of 25 values are averaged to one grid cell.
    if (is.matrix(obj)) {
        obj <- list(x = 1:nrow(obj), y = 1:ncol(obj), z = obj)
    }
    M <- length(obj$x)
    N <- length(obj$y)
    Mi <- trunc(M/Q)
    Ni <- trunc(N/Q)
    # space to hold results
    z <- matrix(NA, nrow = Mi, ncol = N)
    x2 <- rep(NA, Mi)
    y2 <- rep(NA, Ni)
    indQ <- 1:Q
    # sum over block of rows and handle x grid values
    for (j in 1:Mi) {
        x2[j] <- mean(obj$x[indQ + (j - 1) * Q])
        z[j, ] <- colMeans(obj$z[indQ + (j - 1) * Q, ], na.rm = TRUE)
    }
    # sum over blocks of columns  and average y grid values
    for (k in 1:Ni) {
        y2[k] <- mean(obj$y[indQ + (k - 1) * Q])
        z[, k] <- rowMeans(z[, indQ + (k - 1) * Q], na.rm = TRUE)
    }
    return(list(x = x2, y = y2, z = z[1:Mi, 1:Ni], Q = Q))
}
"get.rectangle" <- function() {
    temp <- locator(2, type = "p", pch = "+")
    rect(temp$x[1], temp$y[1], temp$x[2], temp$y[2])
    temp
}
"half.image" <- function(obj) {
    # coerce to list if a matrix
    if (is.matrix(obj)) {
        obj <- list(x = 1:nrow(obj), y = 1:ncol(obj), z = obj)
    }
    M <- length(obj$x)
    N <- length(obj$y)
    M2 <- trunc(M/2)
    N2 <- trunc(N/2)
    z <- matrix(NA, nrow = M2, ncol = N2)
    ix <- (1:M2) * 2
    iy <- (1:N2) * 2
    x2 <- (obj$x[ix - 1] + obj$x[ix])/2
    y2 <- (obj$y[iy - 1] + obj$y[iy])/2
    return(list(x = x2, y = y2, z = (obj$z[ix - 1, iy] + obj$z[ix - 
        1, iy - 1] + obj$z[ix, iy - 1] + obj$z[ix, iy])/4))
}
pushpin <- function(x, y, z, p.out, height = 0.05, 
    col = "black", text = NULL, adj = -0.1, cex = 1, ...) {
    Sxy1 <- trans3d(x, y, z, p.out)
    #trans3d( x,y,z, p.out)-> Sxy2
    Sxy2 <- Sxy1
    hold <- par()$usr
    Sxy2$y <- (hold[4] - hold[3]) * height + Sxy2$y
    segments(Sxy1$x, Sxy1$y, Sxy2$x, Sxy2$y, col = "black")
    points(Sxy2, col = col, pch = 19, cex = cex)
    if (!is.null(text)) {
        text(Sxy2$x, Sxy2$y, label = text, adj = adj, cex = cex, 
            ...)
    }
}
designer.colors <- function(n = 256, col = c("darkgreen", 
    "white", "darkred"), x = seq(0, 1, , length(col))) {
    # distribute the colors at equal spacing if q is NULL
    xg <- seq(0, 1, , n)
    # matrix to hold RGB color values
    y.rgb <- t(col2rgb(col))
    temp <- matrix(NA, ncol = 3, nrow = n)
    # spline interpolation of color values
    for (k in 1:3) {
        hold <- splint(x, y.rgb[, k], xg)
        # fix up to be integer in [0,255]
        hold[hold < 0] <- 0
        hold[hold > 255] <- 255
        temp[, k] <- round(hold)
    }
    # convert back to hex
    rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
}
#boulder.colors<- c('darkred', 'darkorange',
#                   'white', 'darkgreen', 'darkblue')
"two.colors" <- function(n = 256, start = "darkgreen", 
    end = "red", middle = "white") {
    designer.colors(n, c(start, middle, end))
}
