"Krig.find.gcvmin" <-
function(info, lambda.grid, gcv, gcv.fun, tol, verbose = FALSE)
{
	# first take out NA from gcv
	ind <- !is.na(gcv)
	lambda.grid <- lambda.grid[ind]
	gcv <- gcv[ind]
	#
	#find coarse minima
	#
	nstep.cv <- length(lambda.grid)
	il <- order(gcv)[1]
	lambda.gcv <- lambda.grid[il]
	#
	#
	#
	gcv.raw <- min(gcv)
	if(verbose) {
		cat("GCV coarse search:", gcv.raw)
	}
	#
	#
	# do a golden section refined search for minimizing lamdda
	# if the minimum is in interior of the grid search. 
	#
	if((il > 1) & (il < nstep.cv)) {
		#
		# now do the Golden section refinement
		# tolerance for convergence scaled with respect to GCV from the coarse search
		#
		out <- golden.section.search(lambda.grid[il - 1], lambda.grid[
			il], lambda.grid[il + 1], gcv.fun, f.extra = info,
			tol = tol * gcv.raw)
		return(out$x)
	}
	else {
		warning("GCV search gives a minumum at the endpoints of the grid search"
			)
		return(lambda.gcv)
	}
}
