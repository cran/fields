"predict.surface" <-
function (object, grid.list = NA, extrap = FALSE, chull.mask, model = 
NA, 
    nx = 30, ny = 30, ...) 
{
    out<- object # hack S3
    if ((length(grid.list) == 1) | (is.na(grid.list)[1])) {
        if (is.null(out$x)) 
            stop("Need a an X matrix in the\noutput object")
        grid.list <- as.list(rep("c", ncol(out$x)))
        grid.list[[1]] <- "x"
        grid.list[[2]] <- "y"
        if (is.null(dimnames(out$x))) {
            temp <- paste("X", 1:ncol(out$x), sep = "")
        }
        else {
            temp <- dimnames(out$x)[[2]]
        }
        names(grid.list) <- temp
    }

    if (is.null(out$x)) 
        xg <- make.surface.grid(grid.list, nx = nx, ny = ny)
    else xg <- make.surface.grid(grid.list, X = out$x, nx = nx, 
        ny = ny)

# at this point xg is the grid for evaluation

    out2 <- as.surface(xg, predict(out, xg, model = model, ...))

    if (!extrap) {
# columns of variables that are the grid
            ind <- c(attr(xg, "format")[, 1])
        if (missing(chull.mask)) {
            chull.mask <- unique.matrix(out$x[, ind])
        }
# set to NA any point in 2d grid that is outside convex hull        
        out2$z[
        in.poly(xg[,ind], xp = chull.mask,convex.hull=TRUE) == 0
               ] <- NA
    }
    out2
}
