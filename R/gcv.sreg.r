"gcv.sreg" <-
function(out, lambda.grid = NA, cost = 1, nstep.cv = 20, rmse = NA, offset = 0,
	trmin = NA, trmax = NA, verbose = T, tol = 1.0000000000000001e-05,
	find.min = T, method = "GCV")
{
	shat.pure.error <- out$shat.pure.error
	pure.ss <- out$pure.ss
	nt <- 2
	np <- out$np
	N <- out$N
	out$cost <- cost
	out$offset <- offset
	if(is.na(trmin))
		trmin <- 2.0499999999999998
	if(is.na(trmax))
		trmax <- out$np * 0.94999999999999996
	#
	# create a reasonable grid for the GCV search if not supplied
	#
	if(is.na(lambda.grid[1])) {
		l2 <- sreg.df.to.lambda(trmax, out$x, out$wt)
		l1 <- sreg.df.to.lambda(trmin, out$x, out$wt)
		lambda.grid <- exp(seq(log(l2), log(l1),  , nstep.cv))
	}
	#
	# done with finding a good default range for lambda
	#
	#
	# begin grid evaluation over lambda
	#
	if(verbose) {
		cat(l1, l2, fill = T)
	}
	nl <- length(lambda.grid)
	V <- V.model <- V.one <- trA <- MSE <- RSS.model <- rep(NA, nl)
	for(k in 1:nl) {
		temp <- sreg.fit(lambda.grid[k], out)
		RSS.model[k] <- temp$rss
		trA[k] <- temp$trace
		V[k] <- temp$gcv
		V.one[k] <- temp$gcv.one
		V.model[k] <- temp$gcv.model
	}
	RSS <- RSS.model + pure.ss
	#
	# The following three versions all agree if the are unique observations
	# and a full spline basis 
	#
	# V.one the real leave-one-out GCV even if there are replicates etc.
	# V.model This is GCV based on collapsing to rep group means 
	#
	shat <- sqrt(RSS/(N - trA))
	#
	#
	gcv.grid <- cbind(lambda.grid, trA, V, V.one, V.model, shat)
	dimnames(gcv.grid) <- list(NULL, c("lambda", "trA", "GCV", "GCV.one",
		"GCV.model", "shat"))
	#	gcv.grid <- as.data.frame(gcv.grid)
	if(verbose) {
		print(gcv.grid)
	}
	if(!find.min) {
		return(list(gcv.grid = gcv.grid))
	}
	lambda.est <- matrix(ncol = 4, nrow = 5, dimnames = list(c("GCV", 
		"GCV.model", "GCV.one", "RMSE", "pure error"), c("lambda",
		"trA", "GCV", "shat")))
	#
	# refine various GCV estimates of lambda
	# Krig version is generic enough where it can be used unmodified.
	lambda.est[1, 1] <- Krig.find.gcvmin(out, lambda.grid, gcv.grid[, "GCV"
		], sreg.fgcv, tol = tol, verbose = verbose)
	if(verbose) {
		cat(lambda.est[1, 1], fill = T)
	}
	if(!is.na(shat.pure.error)) {
		lambda.est[2, 1] <- Krig.find.gcvmin(out, lambda.grid, gcv.grid[
			, "GCV.model"], sreg.fgcv.model, tol = tol, verbose = 
			verbose)
		if(verbose) {
			cat(lambda.est[2, 1], fill = T)
		}
	}
	lambda.est[3, 1] <- Krig.find.gcvmin(out, lambda.grid, gcv.grid[, 
		"GCV.one"], sreg.fgcv.one, tol = tol, verbose = verbose)
	#
	#
	#    Find lambda to match RMSE if required
	#
	lambda.rmse <- NA
	lambda.pure.error <- NA
	if(!is.na(rmse)) {
		guess <- max(gcv.grid[gcv.grid[, "shat"] < rmse, "lambda"])
		if(verbose) {
			print(rmse)
			print(guess)
		}
		if(!is.na(guess)) {
			lambda.rmse <- find.upcross(sreg.fs2hat, out, 
				upcross.level = rmse^2, guess = guess, tol = 
				tol * rmse^2)
			lambda.est[4, 1] <- lambda.rmse
		}
		else {
			warning("Value of rmse is outside possible range")
		}
	}
	if(!is.na(shat.pure.error)) {
		guess <- max(gcv.grid[gcv.grid[, "shat"] < shat.pure.error,
			"lambda"])
		if(!is.na(guess)) {
			lambda.pure.error <- find.upcross(sreg.fs2hat, out,
				upcross.level = shat.pure.error^2, guess = 
				guess, tol = tol * shat.pure.error^2)
			lambda.est[5, 1] <- lambda.pure.error
		}
		else {
			warning("Value of pure error estimate  \nis outside possible range"
				)
		}
	}
	#
	##
	#
	if(verbose) {
		print(lambda.est)
	}
	for(k in 1:5) {
		lam <- lambda.est[k, 1]
		if(!is.na(lam)) {
			temp <- sreg.fit(lam, out)
			lambda.est[k, 2] <- temp$trace
			if((k == 1) | (k > 3)) {
				lambda.est[k, 3] <- temp$gcv
			}
			if(k == 2) {
				lambda.est[k, 3] <- temp$gcv.model
			}
			if(k == 3) {
				lambda.est[k, 3] <- temp$gcv.one
			}
			lambda.est[k, 4] <- temp$shat
		}
	}
	list(gcv.grid = gcv.grid, lambda.est = lambda.est, lambda.best = 
		lambda.est[method, 1])
}
