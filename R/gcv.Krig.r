"gcv.Krig" <-
function(out, lambda.grid = NA, cost = 1, nstep.cv = 80, rmse = NA, verbose = F,
	tol = 1.0000000000000001e-05, offset = 0, y = NULL, lambda = NA)
{
	nt <- out$nt
	np <- out$np
	N <- out$N
	D <- out$matrices$D
	#
	# If new y's have been passed then update the u vector and the 
	# pure error estimates. Otherwise just use what is in the Krig object
	#
	#
	#
	#
	if(is.null(y)) {
		u <- out$matrices$u
		shat.pure.error <- out$shat.pure.error
		pure.ss <- out$pure.ss
	}
	else {
		# updating part 
		out2 <- Krig.updateY(out, y)
		u <- out2$u
		shat.pure.error <- out2$shat.pure.error
		pure.ss <- out2$pure.ss
	}
	#
	# create a reasonable grid for the GCV search if not supplied
	#
	if(is.na(lambda.grid[1])) {
		l1 <- 1/D[np - nt - 1]
		tr <- np
		#
		########## find upper value of lambda
		#
		for(k in 1:20) {
			tr <- sum(1/(1 + l1 * D))
			if(tr < (nt + 0.050000000000000003))
				break
			l1 <- l1 * 4
		}
		#
		########## find lower lambda
		#
		l2 <- 1/max(D)
		for(k in 1:20) {
			tr <- sum(1/(1 + l2 * D))
			cut <- (N - (cost * (tr - nt - offset) + nt) * (N/
				length(D)))
			if((tr > (np - offset) * 0.94999999999999996) | (cut/
				N <= 0.050000000000000003))
				break
			l2 <- l2/4
		}
		lambda.grid <- exp(seq(log(l2), log(l1),  , nstep.cv))
	}
	#
	# done with finding a good default range for lambda
	#
	#
	# begin grid evaluation over lambda
	#
	nl <- length(lambda.grid)
	nd <- length(D)
	trA <- MSE <- RSS.model <- rep(NA, nl)
	Dl <- rep(NA, nd)
	for(k in 1:nl) {
		Dl <- lambda.grid[k] * D
		RSS.model[k] <- sum(((u * Dl)/(1 + Dl))^2)
		trA[k] <- sum(1/(1 + Dl))
	}
	RSS <- pure.ss + RSS.model
	MSE <- pure.ss/N + RSS.model/length(D)
	MSE.one <- (pure.ss + RSS.model)/N
	den <- 1 - (cost * (trA - nt - offset) + nt)/length(D)
	den.one <- 1 - (cost * (trA - nt - offset) + nt)/N
	#
	# If the denominator is negative then flag this as a bogus case
	# by making the GCV function NA
	#
	#
	# The following three versions all agree if the are unique observations
	#and a 
	# a full basis
	#
	# in the case of replicates this is leave-one-rep_group-out GCV
	V <- ifelse(den > 0, (MSE)/den^2, NA)
	# This is real leave-one-out GCV even if there are replicates etc.
	V.one <- ifelse(den.one > 0, MSE.one/den.one^2, NA)
	#
	# This is GCV based on collapsing to rep group means or in the case of 
	# basis splines prediction relative to the staturated model (lambda=0)
	#
	V.model <- ifelse(den > 0, ((RSS.model/length(D))/den^2), NA)
	#
	#
	shat <- sqrt(RSS/(N - trA))
	#
	#
	# stuff information into a data frame
	#
	gcv.grid <- cbind(lambda.grid, trA, V, V.one, V.model, shat)
	gcv.grid <- as.data.frame(gcv.grid)
	names(gcv.grid) <- c("lambda", "trA", "GCV", "GCV.one", "GCV.model",
		"shat")
	#
	#
	## find global minimum of the GCV function on the grid
	#
	#  create a mini Krig object list with the information needed
	# for further refinements of this estimate and the others
	#
	info <- list(matrices = list(D = D, u = u), N = N, nt = nt, cost = cost,
		pure.ss = pure.ss, shat.pure.error = shat.pure.error, offset = 
		offset)
	#       
	if(verbose) print(info)
	#
	#
	# data frame to hold different estimates for lambda
	#
	lambda.est <- matrix(ncol = 4, nrow = 5, dimnames = list(c("GCV", 
		"GCV.model", "GCV.one", "RMSE", "pure error"), c("lambda",
		"trA", "GCV", "shat")))
	#
	# refine various GCV estimates of lambda
	lambda.est[1, 1] <- Krig.find.gcvmin(info, lambda.grid, gcv.grid$GCV,
		Krig.fgcv, tol = tol, verbose = verbose)
	#
	if(!is.na(shat.pure.error)) {
		temp <- gcv.grid$GCV.model
		lambda.est[2, 1] <- Krig.find.gcvmin(info, lambda.grid, temp,
			Krig.fgcv.model, tol = tol, verbose = verbose)
	}
	lambda.est[3, 1] <- Krig.find.gcvmin(info, lambda.grid, gcv.grid$
		GCV.one, Krig.fgcv.one, tol = tol, verbose = verbose)
	#
	#
	#    Find lambda to match RMSE if required
	#
	lambda.rmse <- NA
	lambda.pure.error <- NA
	if(!is.na(rmse)) {
		guess <- max(gcv.grid$lambda[gcv.grid$shat < rmse])
		if(verbose) {
			print(rmse)
			print(guess)
		}
		if(!is.na(guess)) {
			lambda.rmse <- find.upcross(Krig.fs2hat, info, 
				upcross.level = rmse^2, guess = guess, tol = 
				tol * rmse^2)
			lambda.est[4, 1] <- lambda.rmse
		}
		else {
			warning("Value of rmse is outside possible range")
		}
	}
	#
	##
	#
	if(!is.na(shat.pure.error)) {
		guess <- max(gcv.grid$lambda[gcv.grid$shat < shat.pure.error])
		if(!is.na(guess)) {
			lambda.pure.error <- find.upcross(Krig.fs2hat, info,
				upcross.level = shat.pure.error^2, guess = 
				guess, tol = tol * shat.pure.error^2)
			lambda.est[5, 1] <- lambda.pure.error
		}
		else {
			warning("Value of pure error estimate  is outside possible range"
				)
		}
	}
	#
	#
	# fill in other stuff for each estimate of lambda
	for(k in 1:5) {
		lam <- lambda.est[k, 1]
		if(!is.na(lam)) {
			lambda.est[k, 2] <- Krig.ftrace(lam, D)
			if(k == 1 | k > 3) {
				lambda.est[k, 3] <- Krig.fgcv(lam, info)
			}
			if(k == 2) {
				lambda.est[k, 3] <- Krig.fgcv.model(lam, info)
			}
			if(k == 3) {
				lambda.est[k, 3] <- Krig.fgcv.one(lam, info)
			}
			lambda.est[k, 4] <- sqrt(Krig.fs2hat(lam, info))
		}
	}
	list(gcv.grid = gcv.grid, lambda.est = lambda.est, lambda.best = 
		lambda.est[1, 1])
}
