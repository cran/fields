"summary.gcv.Krig" <-
function(out, lambda, cost = 1, verbose = F, offset = 0, y = NULL)
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
	#  create a mini Krig object list with the information needed
	# for finding statistics
	#
	info <- list(matrices = list(D = D, u = u), N = N, nt = nt, cost = cost,
		pure.ss = pure.ss, shat.pure.error = shat.pure.error, offset = 
		offset)
	#       
	if(verbose) {
		print(info)
	}
	#
	#
	# data frame to hold different estimates for lambda
	#
	lambda.est <- rep(NA, 6)
	names(lambda.est) <- c("lambda", "trA", "GCV", "GCV.one", "GCV.model",
		"shat")
	#
	# fill in stuff for this  lambda
	lambda.est[1] <- lambda
	lambda.est[2] <- Krig.ftrace(lambda, D)
	lambda.est[3] <- Krig.fgcv(lambda, info)
	lambda.est[4] <- Krig.fgcv.one(lambda, info)
	if(!is.na(shat.pure.error)) {
		lambda.est[5] <- Krig.fgcv.model(lambda, info)
	}
	lambda.est[6] <- sqrt(Krig.fs2hat(lambda, info))
	lambda.est
}
