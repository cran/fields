"make.surface.grid" <-
function (grid.list, X, nx = 30, ny = 30, info.list = FALSE, FUN = median) 
{
    what <- rep(NA, 2)
    if (!is.list(grid.list)) 
        stop("Must supply a list to describe grid limits")
    if (!missing(X)) {
        if (data.class(X) != "data.frame") {
            if (is.null(dimnames(X))) {
                names.of.X <- paste("X", 1:ncol(X), sep = "")
            }
            else {
                names.of.X <- dimnames(X)[[2]]
            }
        }
        else names.of.X <- names(X)
        m <- length(names.of.X)
        if (is.null(names(grid.list))) {
            if (length(grid.list) < m) 
                stop(" grid.list must be as long as the number of columns of X!")
            names(grid.list) <- paste("X", 1:length(grid.list), 
                sep = "")
        }
        test <- match(names(grid.list), names.of.X)
        if (!(all(!is.na(test)))) {
            print("names in grid.list")
            print(names(grid.list))
            print("names for columns of X matrix")
            print(names.of.X)
            stop(" some of the\ngrid.list names are not found in the names of the X columns")
        }
        temp <- as.list(rep("c", m))
        names(temp) <- names.of.X
        temp[names(grid.list)] <- grid.list
        for (k in 1:length(temp)) {
            test <- temp[[k]]
            if (length(test) == 1) {
                if (test == "c") 
                  temp[[k]] <- FUN(X[, k])
                if (test == "x") {
                  temp[[k]] <- seq(min(X[, k]), max(X[, k]), 
                    , nx)
                  what[1] <- k
                }
                if (test == "y") {
                  temp[[k]] <- seq(min(X[, k]), max(X[, k]), 
                    , ny)
                  what[2] <- k
                }
            }
        }
        grid.list <- temp
    }
    ind <- unlist(lapply(grid.list, length))
    if (sum(ind > 1) > 2) {
        stop("Only two components can have more than one\nvalue in the grid list")
    }
    nl <- length(grid.list)
    if (is.na(what[1])) {
        what <- (1:nl)[ind > 1]
    }
    x1 <- grid.list[[what[1]]]
    x2 <- grid.list[[what[2]]]
    if (length(x1) == 2) {
        x1 <- seq(min(x1), max(x1), , nx)
    }
    if (length(x2) == 2) {
        x2 <- seq(min(x2), max(x2), , ny)
    }
    nx <- length(x1)
    ny <- length(x2)
    nr <- nx * ny
    if (!info.list) {
        xg <- matrix(NA, ncol = nl, nrow = nr)
        dimnames(xg) <- list(NULL, names(grid.list))
        tempa <- attributes(xg)
        tempa$format <- cbind(what, ind[what])
        tempa$surface.info <- list(x = x1, y = x2, nx = nx, ny = ny, 
            xlab = names(grid.list)[what[1]], ylab = names(grid.list)[what[2]], 
            fixed.variables = grid.list[-what], nvar = nl)
        tempa$grid.list <- grid.list
        attributes(xg) <- tempa
        xg[, what] <- cbind(rep(x1, ny), rep(x2, rep(nx, ny)))
        for (k in 1:nl) {
            if (ind[k] == 1) {
                xg[, k] <- rep(grid.list[[k]], nr)
            }
        }
        class(xg) <- "surface.grid"
        return(xg)
    }
    else {
        return(list(x = x1, y = x2, nx = nx, ny = ny, xlab = names(grid.list)[what[1]], 
            ylab = names(grid.list)[what[2]], fixed.variables = grid.list[-what], 
            grid.list = grid.list, nvar = nl))
    }
}
