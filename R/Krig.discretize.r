"Krig.discretize" <-
function(x, m = 64, n = 64, xg = NULL, yg = NULL, grid = NULL, expand = c(
	1, 1))
{
	#
	# set up discretized grid based on x
	#
	out <- list()
	if(!is.null(grid)) {
		xg <- grid$x
		yg <- grid$y
		m <- length(xg)
		n <- length(yg)
	}
	if(length(expand) == 1)
		expand <- rep(expand, 2)
	if(is.null(xg)) {
		xr <- range(x[, 1])
		deltemp <- (xr[2] - xr[1]) * (expand[1] - 1) * 0.5
		xg <- seq(xr[1] - deltemp, xr[2] + deltemp,  , m)
	}
	else {
		xr <- range(xg)
	}
	if(is.null(yg)) {
		yr <- range(x[, 2])
		deltemp <- (yr[2] - yr[1]) * (expand[2] - 1) * 0.5
		yg <- seq(yr[1] - deltemp, yr[2] + deltemp,  , n)
	}
	else {
		yr <- range(yg)
	}
	del <- 1e-08 + (0.5 * (max(xr[2]) - min(xr[1])))/(m - 1)
	xcut <- seq(xr[1] - del, xr[2] + del,  , (m + 1))
	del <- 1e-08 + (0.5 * (max(yr[2]) - min(yr[1])))/(n - 1)
	ycut <- seq(yr[1] - del, yr[2] + del,  , (n + 1))
	out$m <- length(xg)
	out$n <- length(yg)
	out$index <- cbind(cut(x[, 1], xcut), cut(x[, 2], ycut))
	#	out$cut.grid <- list(x = xcut, y = ycut)	#
	#
	out$grid <- list(x = xg, y = yg)
	out$loc <- cbind(xg[out$index[, 1]], yg[out$index[, 2]])
	out
}
