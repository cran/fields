"Krig.coef" <-
function(out, u = NULL, lambda = NULL, yM = NULL)
{
	# check for C argument.
	C.arg.missing <- is.null(formals(get(out$call.name))$C)
	#
	# the key to computing the estimate are the coeffients c and d
	# or if these need to be recomputed, the vector u
	# for new data, u needs to be found 
	#
	if(is.null(u)) {
		u <- out$matrices$u
	}
	if(is.null(lambda)) {
		lambda <- out$lambda
	}
	if(is.null(yM)) {
		temp.yM <- out$yM
	}
	else {
		temp.yM <- yM
	}
	# use a different lambda or new data so we need to get the # new temp.d and
	#   temp.c  coefficients
	nt <- out$nt
	np <- out$np
	if(out$decomp == "DR") {
		beta <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$
			D)) * u)
		temp.d <- beta[1:nt]
		temp <- c(rep(0, nt), beta[(nt + 1):np])
		#
		#
		# tranform the beta into the parameter associated with the covariance
		# function  basis set 
		#
		temp.c <- c(qr.qy(out$matrices$qr.T, temp))
	}
	if(out$decomp == "WBW") {
		xc <- out$transform$x.center
		xs <- out$transform$x.scale
		knots <- scale(out$knots, xc, xs)
		beta <- out$matrices$G %*% ((1/(1 + lambda * out$matrices$
			D)) * u)
		temp.c <- c(qr.qy(out$matrices$qr.T, c(rep(0, nt), beta[(nt +
			1):np])))
		temp.c <- temp.c * sqrt(out$weightsM)
		# end of if C.arg.missing stmt
		if(C.arg.missing) {
			if(!out$cov.by.name) {
				temp <- temp.yM - lambda * temp.c - out$
					cov.function(knots, knots) %*% temp.c
			}
			# end of if !out$cov.by.name stmt
			if(out$cov.by.name) {
				temp <- temp.yM - lambda * temp.c - do.call(
					out$call.name, c(out$args, list(x1 = 
					knots, x2 = knots))) %*% temp.c
			}
		}
		else {
			if(!out$cov.by.name) {
				temp <- temp.yM - lambda * temp.c - out$
					cov.function(knots, knots, C = temp.c)
			}
			# end of if !out$cov.by.name stmt
			if(out$cov.by.name) {
				temp <- temp.yM - lambda * temp.c - do.call(
					out$call.name, c(out$args, list(x1 = 
					knots, x2 = knots, C = temp.c)))
			}
		}
		# end of if else C.arg.missing stmt
		# multiply through by weights
		temp <- sqrt(out$weightsM) * temp
		temp.d <- qr.coef(out$matrices$qr.T, temp)
	}
	return(list(c = temp.c, d = temp.d))
}
