"predict.Krig" <-
function(out, x = NULL, lambda = NA, df = NA, model = NA, 
	eval.correlation.model = T, y = NULL, verbose = F, gcv = F)
{
	#
	# the key to computing the estimate are the coeffients c and d
	# or if these need to be recomputed, the vector u
	# for new data, u needs to be found 
	#
	if(is.null(y)) {
		temp.u <- out$matrices$u
		temp.yM <- out$yM
		temp.c <- out$c
		temp.d <- out$d
	}
	else {
		out2 <- Krig.updateY(out, y, verbose = verbose)
		temp.u <- out2$u
		temp.yM <- out2$yM
		if(verbose) {
			cat("u's")
			print(temp.u)
		}
	}
	#
	# temp.c and temp.d will be calculated below from the new stuff in out2
	#
	#
	# use data x matrix for prediction if it has been omitted. 
	#
	if(is.null(x)) {
		x <- out$x
	}
	x <- as.matrix(x)
	if(verbose)
		print(x)
	#
	# 
	# sometimes one wants to over ride evaluating as a correlation model
	# check for this flag
	#
	correlation.model <- (out$correlation.model & eval.correlation.model)
	# correlation.model is False unless mean.obj and sd.obj 
	# are specified in the call to Krig
	if(correlation.model) {
		#
		# find mean and sd if this is a correlation model 
		# the final predictions will then be transformed by these quantities
		#
		temp.mean <- predict(out$mean.obj, x)
		temp.sd <- predict(out$sd.obj, x)
	}
	#
	# scale the x values 
	# using information from the output object
	# scaling is (0,1) by default
	#
	if(is.null(out$transform)) {
		xc <- rep(0, ncol(x))
		xs <- rep(1, ncol(x))
	}
	else {
		xc <- out$transform$x.center
		xs <- out$transform$x.scale
	}
	x <- scale(x, xc, xs)
	knots <- scale(out$knots, xc, xs)
	#
	if(verbose) {
		print(x)
	}
	#
	#
	# figure out if coefficients  c and d need to be recomputed
	#	
	find.coef <- (!is.null(y) | !is.na(lambda) | !is.na(df))
	#
	# Now we make various choices for filling in lambda
	#
	# if  degrees of freedom is passed then convert to a lmbda and use
	#this.
	if(!is.na(df)) {
		lambda <- Krig.df.to.lambda(df, out$matrices$D)
	}
	#   Fill in from model component
	if(!is.na(model)) {
		lambda <- model[1]
	}
	#
	# if lambda is actually passed and no GCV then use this value for lambda
	#
	if(is.na(lambda) & !gcv) lambda <- out$lambda
	#
	# if GCV then do it and find lambda
	# here we are assuming that new data is supplied.  
	#
	if(gcv) {
		lambda <- gcv.Krig(out, cost = out$cost, offset = out$offset,
			y = y)$lambda.best
	}
	if(find.coef) {
		out3 <- Krig.coef(out, u = temp.u, lambda = lambda, yM = 
			temp.yM)
		temp.d <- out3$d
		temp.c <- out3$c
		if(verbose) {
			cat(" d coefs")
			print(temp.d)
			cat("c coefs")
			print(temp.c)
		}
	}
	#
	# At this point we have  the right coefficients (temp.d and temp.c)
	# to compute the spline.  
	#
	#
	# decide whether to use the fast multiplication routines for the
	#covariance function
	#
	#
	#
	# check whether the covarinace function has the argument C in its call.
	# If so take advantage of it. 
	#
	C.arg.missing <- is.null(formals(get(out$call.name))$C)
	if(C.arg.missing) {
		if(!out$cov.by.name) {
			temp <- c(out$make.tmatrix(x, out$m) %*% temp.d + out$
				cov.function(x, knots) %*% temp.c)
		}
		# end of if !out$cov.by.name
		if(out$cov.by.name) {
			temp <- c(out$make.tmatrix(x, out$m) %*% temp.d + 
				do.call(out$call.name, c(out$args, list(x1 = x,
				x2 = knots))) %*% temp.c)
		}
	}
	else {
		if(!out$cov.by.name) {
			temp <- c(out$make.tmatrix(x, out$m) %*% temp.d + out$
				cov.function(x, knots, C = temp.c))
		}
		# end of if !out$cov.by.name stmt
		if(out$cov.by.name) {
			temp <- c(out$make.tmatrix(x, out$m) %*% temp.d + 
				do.call(out$call.name, c(out$args, list(x1 = x,
				x2 = knots, C = temp.c))))
		}
	}
	#
	# if correlation model do the transformation
	#
	if(correlation.model) temp <- (temp * temp.sd + temp.mean)
	#
	# return lambda also if GCV was done
	#
	if(!gcv) {
		return(temp)
	}
	if(gcv) {
		return(list(predicted = temp, lambda = lambda))
	}
}
