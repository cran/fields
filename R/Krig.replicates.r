"Krig.replicates" <-
function(out, verbose = FALSE)
{
	rep.info <- cat.matrix(out$x)
	if(verbose)
		print(rep.info)
	uniquerows <- !duplicated(rep.info)
	if(sum(uniquerows) == out$N) {
		shat.rep <- NA
		shat.pure.error <- NA
		pure.ss <- 0
		yM <- out$y
		weightsM <- out$weights
		xM <- as.matrix(out$x[uniquerows,  ])
	}
	else {
		rep.info.aov <- fast.1way(rep.info, out$y, out$weights)
		shat.pure.error <- sqrt(rep.info.aov$MSE)
		shat.rep <- shat.pure.error
		yM <- rep.info.aov$means
		weightsM <- rep.info.aov$w.means
		xM <- as.matrix(out$x[uniquerows,  ])
		pure.ss <- rep.info.aov$SSE
		if(verbose)
			print(rep.info.aov)
	}
	return(yM, xM, weightsM, uniquerows, shat.rep, shat.pure.error, pure.ss,
		rep.info)
}
