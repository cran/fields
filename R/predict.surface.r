"predict.surface" <-
function(out, grid.list = NA, extrap = F, chull.mask, model = NA, ...)
{
	if((length(grid.list) == 1) | (is.na(grid.list)[1])) {
		if(is.null(out$x))
			stop("Need a an X matrix in the\noutput object")
		grid.list <- as.list(rep("c", ncol(out$x)))
		grid.list[[1]] <- "x"
		grid.list[[2]] <- "y"
		# add some columns names to the X matrix if there are none
		if(is.null(dimnames(out$x))) {
			temp <- paste("X", 1:ncol(out$x), sep = "")
		}
		else {
			temp <- dimnames(out$x)[[2]]
		}
		names(grid.list) <- temp
	}
	if(is.null(out$x))
		xg <- make.surface.grid(grid.list)
	else xg <- make.surface.grid(grid.list, X = out$x)
	out2 <- as.surface(xg, predict(out, xg, model = model, ...))
	if(!extrap) {
		if(missing(chull.mask)) {
			ind <- c(attr(xg, "format")[, 1])
			chull.mask <- unique.matrix(out$x[, ind])
		}
		#
		# find grid points that are outside convex hull 
		# of 2-d projection of the observations
		# or the set of points passed in the call 
		# set these grid points to NA in the matrix of gridded values
		#               
		chull.mask <- chull.mask[chull(chull.mask),  ]
		out2$z[in.poly(xg, xp = chull.mask) == 0] <- NA
	}
	out2
}
