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
		X <- cbind(make.tmatrix(xM, out$m), qr.yq2(out$matrices$qr.T,
			do.call(out$call.name, c(out$args, list(x1 = xM, x2 = 
			knots)))))
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
			temp <- cbind(make.tmatrix(x0, out$m), matrix(qr.yq2(
				out$matrices$qr.T, do.call(out$call.name, c(
				out$args, list(x1 = x0, x2 = knots)))), nrow = 
				nrow(x0))) %*% temp
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
			tempd <- diag(1, np) - lambda * tempc - do.call(out$
				call.name, c(out$args, list(x1 = knots, x2 = 
				knots))) %*% tempc
		tempd <- diag(sqrt(out$weightsM)) %*% tempd
		tempd <- qr.coef(out$matrices$qr.T, tempd)
		# hat matrix is finding d and c coefficients and 
		# evaluting at basis functions.
			temp <- make.tmatrix(x0, out$m) %*% tempd + do.call(
				out$call.name, c(out$args, list(x1 = x0, x2 = 
				knots))) %*% tempc
	}
	#
	#
	return(temp)
}
