"Krig" <-
function(x, Y, cov.function = "exp.cov", lambda = NA, df = NA, cost = 1., knots,
	weights = rep(1., length(Y)), m = 2., return.matrices = TRUE, nstep.cv = 
	80., scale.type = "user", x.center = rep(0., ncol(x)), x.scale = rep(
	1., ncol(x)), rho = NA, sigma2 = NA, method = "GCV", decomp = "DR",
	verbose = FALSE, cond.number = 10.^8., mean.obj = NULL, sd.obj = NULL,
	yname = NULL, return.X = TRUE, null.function = make.tmatrix, offset = 0.,
	outputcall = NULL, cov.args = NULL, na.rm=FALSE, ...)
{
	# Here is the ouput list for this function
	out <- list()
	class(out) <- c("Krig")

# save the covariance function
	if(!is.character(cov.function)) {
		if(is.function(cov.function))
			cov.function <- as.character(substitute(cov.function))
	}
		out$call.name <- cov.function

# save the call 	
        if(is.null(outputcall)) {
		out$call <- match.call()
	}
	else {
		out$call <- outputcall
	}
#
#
	if(sum(is.na(Y)) != 0.&!na.rm) {
	stop("Need to remove missing values or set na.rm equal to TRUE")
	}

# name of Y data
	if(is.null(yname)) {
             out$yname <- deparse(substitute(Y)) }
        else {
                 out$yname <-yname
        }   
#
# setup up other components of output object.     

	out$decomp <- decomp
	#
	out$make.tmatrix <- null.function
	###
	out$offset <- offset
	out$cost <- cost
		# -- covariance  parameters are passed separately into
		# -- the out list.  They will be stored in out$args
		if(!is.null(cov.args)) out$args <- cov.args else out$args <-
				list(...)
		if(verbose) {
			cat(" Local cov function arguments ", fill = TRUE)
			print(out$args)
			cat(" covariance function used is : ", fill = TRUE)
			print(out$call.name)
		}
	#
	# make sure that the columns of x have names
	#
	# add row and col names if they are missing.
	# 
	x <- as.matrix(x)
	if(length(dimnames(x)) != 2.) {
		dimnames(x) <- list(format(1.:nrow(x)), paste("X", format(
			1.:ncol(x)), sep = ""))
	}
	if(length(dimnames(x)[[1.]]) == 0.) {
		dimnames(x)[[1.]] <- format(1.:nrow(x))
	}
	#
	#
	#
	if(length(dimnames(x)[[2.]]) == 0.) {
		dimnames(x)[[2.]] <- paste("X", format(1.:ncol(x)), sep = "")
	}
	Y <- c(Y)

        #
        # deal with NA's 
        #
        if( na.rm){
        ind<- is.na( Y)
        if( any(ind)){
        Y<- Y[!ind]
        x<- x[!ind,]
        warning("NA's have been removed from Y and corresponding X rows removed.")
        } 
        }

#
# OK here are the real data sets
        N<- length(Y)
	out$N <- N
	out$y <- Y
	out$x <- x
#
	if(!is.null(sd.obj) & !is.null(mean.obj)) {
		correlation.model <- TRUE
	}
	else correlation.model <- FALSE
	out$correlation.model <- correlation.model
	#
	#
	#
	if(correlation.model) {
		Yraw <- Y
		out$mean.obj <- mean.obj
		out$sd.obj <- sd.obj
		Y <- (Y - predict(mean.obj, x))/predict(sd.obj, x)
		if(verbose)
			print(Y)
	}
	#
	## find duplicate rows of the X matrix 
	## first find integer tags to indicate replications
	out$weights <- weights
	rep.info <- cat.matrix(x)
	out$rep.info <- rep.info
	if(verbose) {
		cat("rep.info", fill = TRUE)
		print(rep.info)
	}
	if(max(rep.info) == N) {
		shat.rep <- NA
		shat.pure.error <- NA
		out$pure.ss <- 0.
		YM <- Y
		weightsM <- weights
		xM <- as.matrix(x[!duplicated(rep.info),  ])
	}
	else {
		##
		## do a simple 1-way ANOVA to get the replication error
		##
		rep.info.aov <- fast.1way(rep.info, Y, weights)
		shat.pure.error <- sqrt(rep.info.aov$MSE)
		shat.rep <- shat.pure.error
		YM <- rep.info.aov$means
		weightsM <- rep.info.aov$w.means
		xM <- as.matrix(x[!duplicated(rep.info),  ])
		out$pure.ss <- rep.info.aov$SSE
		if(verbose) {
			cat(" rep info", fill = TRUE)
			print(rep.info.aov)
		}
	}
	out$yM <- YM
	out$xM <- xM
	out$weightsM <- weightsM
	out$shat.rep <- shat.rep
	out$shat.pure.error <- shat.pure.error
	if(missing(knots)) {
		knots <- xM
		mle.calc <- TRUE
		knot.model <- FALSE
	}
	else {
		mle.calc <- FALSE
		knot.model <- TRUE
	}
	out$knot.model <- knot.model
	out$mle.calc <- mle.calc
	knots <- as.matrix(knots)
	##
	#
	## scale x and knots 
	out$knots <- knots
	xM <- transformx(xM, scale.type, x.center, x.scale)
	transform <- attributes(xM)
	if(verbose) {
		cat("transform", fill = TRUE)
		print(transform)
	}
	knots <- scale(knots, center = transform$x.center, scale = transform$
		x.scale)
	if(verbose) {
		cat("knots in transformed scale", fill = TRUE)
		print(knots)
	}
	##
	##
	##  use value of lambda implied by rho and simga2 if these are passed
	##
	out$transform <- transform
	if(!is.na(lambda) | !is.na(df))
		method <- "user"
	if(!is.na(rho) & !is.na(sigma2)) {
		lambda <- sigma2/rho
		method <- "user"
	}
	just.solve <- (lambda[1.] == 0.)
	#
	#
	if(is.na(just.solve)) just.solve <- FALSE
	#
	if(verbose) cat("lambda", lambda, fill = TRUE)
	# find the QR decopmposition of T matrix  that spans null space with
	# respect to the knots 
	# the big X matrix
	#
	d <- ncol(xM)
	if(decomp == "DR") {
		#
		qr.T <- qr(out$make.tmatrix(knots, m))
		if(verbose) {
			print(qr.T)
		}
			tempM <- qr.yq2(qr.T, do.call(out$call.name, c(out$args, list(x1 = knots, x2 = knots))))
		tempM <- qr.q2ty(qr.T, tempM)
		if(verbose) {
			print(dim(tempM))
		}
	}
	if(decomp == "WBW") {
		#
		#   construct the covariance matrix 
		#functions and Qr decomposition of T
		#
		#
		#
		qr.T <- qr(sqrt(weightsM) * out$make.tmatrix(knots, m))
		if(verbose) {
			print(qr.T)
		}
			tempM <- sqrt(weightsM) * t(sqrt(weightsM) * t(do.call(
				out$call.name, c(out$args, list(x1 = knots,
				x2 = knots)))))
		tempM <- qr.yq2(qr.T, tempM)
		tempM <- qr.q2ty(qr.T, tempM)
	}
	np <- nrow(knots)
	# number of para. in NULL space
	nt <- (qr.T$rank)
	out$np <- np
	out$nt <- nt
	if(verbose)
		cat("np, nt", np, nt, fill = TRUE)
	# if lambda = 0 then just solve the system 
	if(verbose) print(knots)
	if(just.solve) {
			beta <- qr.coef(qr(cbind(out$make.tmatrix(xM, m), 
				qr.yq2(qr.T, do.call(out$call.name, c(out$
				args, list(x1 = xM, x2 = knots)))))), YM)
	}
	else {
		#
		#   do all the heavy decompositions if lambda is not = 0
		#   or if it is omitted
		#
		#
		####
		####
		#### Block for Full Demmler Reinsch decomposition
		#####
		#####
		##### end DR decomposition block 
		#####
		####
		#### begin WBW decomposition block
		####
		if(decomp == "DR") {
			## make up penalty matrix
			#
			if(verbose) cat("Type of decomposition", decomp, fill
					 = TRUE)
			H <- matrix(0., ncol = np, nrow = np)
			#
			#
			# svd of big X matrix in preparation for finding inverse square root
			#
			H[(nt + 1.):np, (nt + 1.):np] <- tempM
				X <- cbind(out$make.tmatrix(xM, m), qr.yq2(
					qr.T, do.call(out$call.name, c(out$
					args, list(x1 = xM, x2 = knots)))))
			if(verbose) {
				print(weightsM)
				cat("first 3 rows of X", fill = TRUE)
				print(X[1.:3.,  ])
			}
			temp <- svd(sqrt(weightsM) * X)[c("v", "d")]
			cond.matrix <- max(temp$d)/min(temp$d)
			#
			# inverse symetric square root of X^T W  
			#
			if(cond.matrix > cond.number) {
				print(paste("condition number is ", cond.matrix,
					sep = ""))
				stop("Covariance matrix is close\nto\nsingular"
					)
			}
			#   eigenvalue eigenvector decomposition of BHB
			#
			B <- temp$v %*% diag(1./(temp$d)) %*% t(temp$v)
			temp <- svd(B %*% H %*% B)
			U <- temp$u
			#
			D <- temp$d
			#   We know that H has atleast nt zero singular values ( see how H is
			#   filled)
			#   So make these identically zero.
			#   the singular values are returned from largest to smallest.
			#
			if(verbose) {
				cat("singular values:", fill = TRUE)
				print(D)
			}
			D[(1.:nt) + (np - nt)] <- 0.
			#
			#   with these these decompositions it now follows that 
			#     b= B*U( I + lambda*D)^(-1) U^T * B * X^T*Y
			#      = G*( I + lambda*D)^(-1) G^T* X^T*Y
			#	
			# Now tranform  Y based on this last equation
			#
			G <- B %*% U
			u <- t(G) %*% t(X) %*% (weightsM * YM)
			#
			#
			#   So now we have   
			#
			#    b= G*( I + lambda*D)^(-1)*u 
			#   Note how in this form we can rapidly solve for b for any lambda
			#
			# save matrix decopositions in out list
			#
			# find the pure error sums of sqaures.  
			#
			if(verbose) {
				cat("DR: u", fill = TRUE)
				print(u)
			}
			out$pure.ss <- sum(weightsM * (YM - X %*% G %*% u)^
				2.) + out$pure.ss
			if(verbose) {
				cat("pure.ss", fill = TRUE)
				print(out$pure.ss)
			}
			out$matrices <- list(B = B, U = U, u = u, D = D, G = G,
				qr.T = qr.T, X = X)
		}
		#####
		##### end WBW block
		#####
		###
		### Selection of lambda
		### this is done even if the smoothing parameter or degress of freedom is
		### specified. 
		###
		if(decomp == "WBW") {
			#### decomposition of Q2TKQ2
			temp <- svd(tempM)[c("d", "v")]
			D <- c(rep(0., nt), 1./temp$d)
			if(verbose) {
				cat("singular values:", fill = TRUE)
				print(D)
			}
			G <- matrix(0., ncol = np, nrow = np)
			G[(nt + 1.):np, (nt + 1.):np] <- temp$v
			G <- G * matrix(D, ncol = np, nrow = np, byrow = TRUE)
			u <- c(rep(0., nt), t(temp$v) %*% qr.q2ty(qr.T, sqrt(
				weightsM) * YM))
			#
			# pure error in this case is found for the 1way ANOVA 
			#
			# test for replicates
			#
			if(verbose) {
				cat("WBW: u", fill = TRUE)
				print(u)
			}
			if(verbose) {
				cat("WBW: pure.ss", fill = TRUE)
				print(out$pure.ss)
			}
			out$matrices <- list(u = u, D = D, G = G, qr.T = qr.T,
				decomp = decomp, V = temp$v)
		}
		if(verbose) {
			cat("call to gcv minimization", fill = TRUE)
		}
		gcv.out <- gcv.Krig(out, nstep.cv = nstep.cv, verbose = verbose,
			cost = out$cost, offset = out$offset, lambda = lambda)
		gcv.grid <- gcv.out$gcv.grid
		out$gcv.grid <- gcv.grid
		out$lambda.est <- gcv.out$lambda.est
		if(verbose) {
			print(out$gcv.grid)
		}
	}
	out$method <- method
	if(method == "user") {
		#
		# if degrees of freedom is specified use the this implied value for lambda
		#
		if(!is.na(df)) lambda <- Krig.df.to.lambda(df, D)
		lambda.best <- lambda
		out$lambda.est.user <- summary.gcv.Krig(out, lambda.best, 
			offset = out$offset, cost = out$cost)
	}
	else {
		lambda.best <- gcv.out$lambda.est[method, 1.]
	}
	# add in null space parameters if this WBW decomp
	#
	#
	beta <- G %*% ((1./(1. + lambda.best * D)) * u)
	out$m <- m
	if(!just.solve) {
		out$eff.df <- sum(1./(1. + lambda.best * D))
		out$trace <- out$eff.df
		if(verbose) {
			cat("trace of A", fill = TRUE)
			print(out$trace)
		}
	}
	else {
		out$eff.df <- out$np
	}
	if(just.solve)
		out$lambda <- lambda
	else out$lambda <- lambda.best
	##
	#
	# tranform the beta into the parameter associated with the covariance
	# function
	# basis set. 
	#  into the c parameter vector. 
	#
	out$beta <- beta
	hold <- Krig.coef(out)
	out$c <- hold$c
	# find predicted values and residuals. 
	#
	out$d <- hold$d
	if(verbose) {
		cat(names(out), fill = TRUE)
	}
	#  need eval.correlation.model = F in order to get predicted's in the
	# standardized scale of Y
	# 
	out$fitted.values <- predict.Krig(out, x, eval.correlation.model = FALSE)
	out$residuals <- Y - out$fitted.values
	#
	# funny conversions are in case nt is equal to 1 and X is just a vector
	#
	if(verbose) {
		cat("resdiuals", out$residuals, fill = TRUE)
	}
	#
	#
	#
	out$fitted.values.null <- as.matrix(out$make.tmatrix(x, m)) %*% out$
		d
	#
	out$just.solve <- just.solve
	# fill in the linear parameters of the covariance function in 
	# the output object
	#
	# the next formula is pretty strange. It follows from solving the
	# system of equations for the basis coefficients. 
	#       
	out$shat.GCV <- sqrt(sum(out$weights * out$residuals^2.)/(length(Y) -
		out$trace))
	if(mle.calc) {
		out$rhohat <- sum(out$c * out$yM)/(N - nt)
		if(is.na(rho)) {
			out$rho <- out$rhohat
		}
		else {
			out$rho <- rho
		}
		if(is.na(sigma2))
			sigma2 <- out$rho * out$lambda
		#
		out$sigma2 <- sigma2
		out$shat.MLE <- sqrt(out$rhohat * out$lambda)
		out$best.model <- c(out$lambda, out$sigma2, out$rho)
	}
	## wipe out big matrices if they are not to be returned
	##
	##
	if(!mle.calc) {
		out$sigma2 <- out$shat.GCV^2.
		out$rho <- out$sigma2/out$lambda
		out$rhohat <- NA
		out$shat.MLE <- NA
		out$best.model <- c(out$lambda, out$sigma2)
		out$warning <- 
			"Maximum likelihood estimates not found with knots \n"
		print(paste(out$warning, " \n "))
	}
	if(!return.matrices) {
		out$xM <- NULL
		out$YM <- NULL
		out$x <- NULL
		out$y <- NULL
		out$matrices <- NULL
		out$weights <- NULL
	}
	##
	##
	if(!return.X) out$matrices$X <- NULL
	out
}
