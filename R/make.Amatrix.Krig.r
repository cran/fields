"make.Amatrix.Krig" <-
function(out, x0 = out$x, lambda)
{
	if(missing(lambda)) {
		lambda <- out$lambda
	}
	xc <- out$transform$x.center
	xs <- out$transform$x.scale
	xM <- scale(out$xM, xc, xs)
	knots <- scale(out$knots, xc, xs)
	x0 <- scale(x0, xc, xs)
	if(!out$cov.by.name) {
		X <- cbind(make.tmatrix(xM, out$m), qr.yq2(out$matrices$qr.T,
			out$cov.function(xM, knots)))
	}
	# end of if !out$cov.by.name stmt
	if(out$cov.by.name) {
		X <- cbind(make.tmatrix(xM, out$m), qr.yq2(out$matrices$qr.T,
			do.call(out$call.name, c(out$args, list(x1 = xM, x2 = 
			knots)))))
	}
	# end of if out$cov.by.name stmt
	if(out$decomp == "DR") {
		temp <- (out$matrices$G) %*% diag(1/(1 + lambda * out$matrices$
			D))
		temp <- temp %*% t(out$matrices$G) %*% t(X)
		temp <- temp %*% diag(out$weightsM)
		#
		#
		# At this point temp maps YM 
		# to the beta vector of coefficients
		#
		if(!out$cov.by.name) {
			temp <- cbind(make.tmatrix(x0, out$m), matrix(qr.yq2(
				out$matrices$qr.T, out$cov.function(x0, knots)),
				nrow = nrow(x0))) %*% temp
		}
		# end of if !out$cov.by.name stmt
		if(out$cov.by.name) {
			temp <- cbind(make.tmatrix(x0, out$m), matrix(qr.yq2(
				out$matrices$qr.T, do.call(out$call.name, c(
				out$args, list(x1 = x0, x2 = knots)))), nrow = 
				nrow(x0))) %*% temp
		}
	}
	if(out$decomp == "WBW") {
		nt <- out$nt
		np <- out$np
		temp <- matrix(0, out$np, out$np)
		#
		# matrix that gives the u vector
		temp[(nt + 1):np,  ] <- t(out$matrices$V) %*% qr.q2ty(out$
			matrices$qr.T, diag(sqrt(out$weightsM)))
		#
		# matrix that gives the beta vector
		temp <- out$matrices$G %*% diag((1/(1 + lambda * out$matrices$
			D))) %*% temp
		#
		#		matrix that gives the c vector
		tempc <- diag(sqrt(out$weightsM)) %*% qr.qy(out$matrices$qr.T,
			temp)
		#matrix to get the d vector
		if(!out$cov.by.name) {
			tempd <- diag(1, np) - lambda * tempc - out$
				cov.function(knots, knots) %*% tempc
		}
		# end of if !out$cov.by.name stmt
		if(out$cov.by.name) {
			tempd <- diag(1, np) - lambda * tempc - do.call(out$
				call.name, c(out$args, list(x1 = knots, x2 = 
				knots))) %*% tempc
		}
		# end of if out$cov.by.name stmt
		tempd <- diag(sqrt(out$weightsM)) %*% tempd
		tempd <- qr.coef(out$matrices$qr.T, tempd)
		# hat matrix is finding d and c coefficients and 
		# evaluting at basis functions.
		if(!out$cov.by.name) {
			temp <- make.tmatrix(x0, out$m) %*% tempd + out$
				cov.function(x0, knots) %*% tempc
		}
		# end of if !out$cov.by.name stmt
		if(out$cov.by.name) {
			temp <- make.tmatrix(x0, out$m) %*% tempd + do.call(
				out$call.name, c(out$args, list(x1 = x0, x2 = 
				knots))) %*% tempc
		}
	}
	#
	#
	return(temp)
}
