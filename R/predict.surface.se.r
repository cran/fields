"predict.surface.se" <-
function (object, grid.list = NA, extrap = FALSE, chull.mask, nx=80, ny=80,
    ...) 
{
    out <- object
    if ((length(grid.list) == 1) | (is.na(grid.list)[1])) {
        grid.list <- as.list(rep("c", ncol(out$x)))
        grid.list[[1]] <- "x"
        grid.list[[2]] <- "y"
        temp <- dimnames(out$x)[[2]]
        if (!(temp[1] == "")) 
            names(grid.list) <- temp
    }
    if (is.null(out$x)) 
        xg <- make.surface.grid(grid.list, nx=nx, ny=ny)
    else xg <- make.surface.grid(grid.list, X = out$x, nx=nx, ny=ny)
    out2 <- as.surface(xg, predict.se(out, xg, ...))
    if (!extrap) {
        if (missing(chull.mask)) {
            ind <- c(attr(xg, "format")[, 1])
            chull.mask <- out$x[, ind]
        }
        chull.mask <- chull.mask[chull(chull.mask), ]
        out2$z[in.poly(xg, xp = chull.mask) == 0] <- NA
    }
    out2
}

