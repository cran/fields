"Krig.newy" <-
function(out, newy, verbose = F)
{
	#
	# calculates the collapsed y (weighted) mean vector based on the
	# X matrix and weights from the out object. 
	#
	if(length(newy) != out$N) {
		stop(" the new Y vector is the wrong length!")
	}
	#
	# If there are no replicated obs. then return the full vector
	# pure error ss is zero 
	#
	if(length(unique(out$rep.info)) == out$N) {
		shat.rep <- NA
		shat.pure.error <- NA
		pure.ss <- 0
		yM <- ynew
		return(YM = ynew, shat.rep = NA, shat.pure.error = NA, pure.ss =
			0)
	}
	else {
		#
		# calculate means by pooling Replicated obseravations but use the
		# the right weighting. 
		#
		rep.info.aov <- fast.1way(out$rep.info, ynew, out$weights)[
			c("means", "MSE", "SSE")]
		shat.pure.error <- sqrt(rep.info.aov$MSE)
		shat.rep <- shat.pure.error
		return(yM = rep.info.aov$means, shat.rep = shat.rep, 
			shat.pure.error = shat.pure.error, pure.ss = 
			rep.info.aov$SSE)
	}
}
