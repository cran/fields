"predict.se.Krig" <-
function(out, x, cov.function, rho, sigma2, weights = NULL, cov = F, stationary
	 = T, fixed.mean = T, ...)
{
	# NOTE: to follow along 
	# if lambda is an estimate of the random field at a point
	# and lambda_hat is an estimate based on linear combination of the data
	#
	# MSE( lambda - lambda_hat) = VAR( lambda) - 2* COV( lambda, lambda_hat)
	#      + VAR( lambda_hat)
	#
	# Perhaps the most confusing aspect of this code is that the MSE
	# are found simultanouesly for all the points in x.
	# This is to avoid an explicit for loop 
	#
	#
	if(missing(x)) x <- out$x
	x <- as.matrix(x)
	xraw <- x
	#### transform x to the scale used to fit model
	if(is.null(out$transform)) {
		xc <- rep(0, ncol(x))
		xs <- rep(1, ncol(x))
	}
	else {
		xc <- out$transform$x.center
		xs <- out$transform$x.scale
	}
	x <- scale(x, xc, xs)
	xM <- scale(out$xM, xc, xs)
	##
	########## if covariance function and parameters are missing
	########## extract them from the Krig object
	##
	## figure out what covariance function to use. 
	if(!out$cov.by.name) {
		if(missing(cov.function)) {
			fun <- out$cov.function
		}
		else {
			fun <- cov.function
		}
	}
	# end of if !out$cov.by.name
	if(missing(sigma2)) {
		sigma2 <- out$sigma2
	}
	#
	# OK now fix the right value for sigma
	#
	if(missing(rho)) {
		rho <- out$rho
	}
	#
	# set switch if the covariance is stationary
	if(!out$cov.by.name) {
		if(!is.null(fun$marginal))
			stationary <- F
	}
	# end of if !out$cov.by.name stmt
	if(out$cov.by.name) {
		if(!is.null(formals(get(out$call.name))$marginal)) {
			stationary <- F
		}
	}
	# end of if out$cov.by.name stmt
	#
	lambda <- sigma2/rho
	if(lambda != 0) {
		if((lambda - out$lambda)/lambda > 1e-08) {
			warning("lambda value used is different from the one in the Krig object"
				)
			cat(lambda, out$lambda, fill = T)
		}
	}
	nx <- nrow(xM)
	#
	#       figure out if this is a correlation model and make adustments
	#
	temp.sd <- 1
	var.mean <- 0
	if(out$correlation.model) {
		# predict standard deviations at xraw values
		temp.sd <- c(predict(out$sd.obj, xraw))
		if(!fixed.mean) {
			# for the very end ... the variance in estimating the mean function. 
			var.mean <- predict.se(out$mean.obj, xraw)^2
		}
	}
	#
	# here is where the real fun begins
	#
	# wght.vec are the linear combinations of the data ( yM the averages
	# if there are replicates) that give the 
	# correpsonding estimates of the function at the points x
	wght.vec <- t(make.Amatrix(out, xraw, lambda))
	#
	# modify weight.vec to give the weighted estimated 
	#	
	# Cy is the observed covariance matrix of the data vector
	# xM should be in transformed scale and are the unique x values from
	# the data
	if(!out$cov.by.name) {
		Cy <- rho * fun(xM, xM, ...) + sigma2 * diag(1/out$weightsM)
	}
	# end of if !out$cov.by.name stmt
	if(out$cov.by.name) {
		Cy <- rho * do.call(out$call.name, c(out$args, list(x1 = xM,
			x2 = xM))) + sigma2 * diag(1/out$weightsM)
	}
	# end of if out$cov.by.name stmt
	#
	# see if the C argument is given, if so utilize it.
	C.arg.missing <- is.null(formals(get(out$call.name))$C)
	if(cov) {
		# this block find the whole covariance matrix 
		if(!out$cov.by.name) {
			temp <- rho * fun(x, xM, ...) %*% wght.vec
			temp <- (rho * fun(x, x, ...) - temp - t(temp))
		}
		# end of if !out$cov.by.name stmt
		if(out$cov.by.name) {
			if(C.arg.missing) {
				temp <- rho * do.call(out$call.name, c(out$
					args, list(x1 = x, x2 = xM))) %*% 
					wght.vec
			}
			else {
				temp <- rho * do.call(out$call.name, c(out$
					args, list(x1 = x, x2 = xM, C = 
					wght.vec)))
			}
			# end of if else C.arg.missing stmt
			temp <- rho * do.call(out$call.name, c(out$args, list(
				x1 = x, x2 = x))) - temp - t(temp)
		}
		# end of if out$cov.by.name stmt
		temp <- temp + t(wght.vec) %*% Cy %*% wght.vec
		return(temp.sd * t(t(temp * temp.sd)))
	}
	#
	# this is the second term in the MSE
	if(!out$cov.by.name) {
		temp1 <- rho * c(t(wght.vec * fun(xM, x, ...)) %*% rep(1, nx))
	}
	# end of if !out$cov.by.name stmt
	if(out$cov.by.name) {
		temp1 <- rho * c(t(wght.vec * do.call(out$call.name, c(out$
			args, list(x1 = xM, x2 = x)))) %*% rep(1, nx))
	}
	# end of if out$cov.by.name stmt
	#
	# this is the third term in the MSE 
	temp2 <- c(t(wght.vec * (Cy %*% wght.vec)) %*% rep(1, nx))
	#
	# VAR of the first term is handled differently depending on 
	# if this is the MSE at a point or a weighted combination
	#
	if(stationary) {
		x0 <- matrix(0, ncol = ncol(x), nrow = 1)
		if(!out$cov.by.name) {
			temp <- (rho * fun(x0, x0, ...) - 2 * temp1 + temp2)
		}
		# end of if !out$cov.by.name stmt
		if(out$cov.by.name) {
			temp <- rho * do.call(out$call.name, c(out$args, list(
				x1 = x0, x2 = x0))) - 2 * temp1 + temp2
		}
		# end of if out$cov.by.name stmt
		return(sqrt(temp * temp.sd^2 + var.mean))
	}
	#
	#
	# if covariance is not stationary then loop through each point to get
	# the variance of field at that point. 
	#
	#
	if(!out$cov.by.name) {
		if(is.null(fun$marginal)) {
			temp <- rep(0, nrow(x))
			for(k in 1:nrow(x)) {
				x0 <- matrix(x[k,  ], nrow = 1)
				temp[k] <- rho * fun(x0, x0, ...) - 2 * temp1[
					k] + temp2[k]
			}
		}
		else {
			#
			## marginal variances available by single call
			#
			temp <- rho * fun(x, marginal = T, ...) - 2 * temp1 +
				temp2
		}
	}
	# end of if !out$cov.by.name stmt
	if(out$cov.by.name) {
		if(is.null(formals(get(out$call.name))$marginal)) {
			temp <- rep(0, nrow(x))
			for(k in 1:nrow(x)) {
				x0 <- matrix(x[k,  ], nrow = 1)
				temp[k] <- rho * do.call(out$call.name, c(
					out$args, list(x1 = x0, x2 = x0))) -
					2 * temp1[k] + temp2[k]
			}
		}
		else {
			#
			## marginal variances available by single call
			#
			temp <- rho * do.call(out$call.name, c(out$args, list(
				x, marginal = T))) - 2 * temp1 + temp2
		}
	}
	# end of if out$cov.by.name stmt
	return(sqrt((temp.sd^2) * temp + var.mean))
}
