# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"gcv.sreg" <-
function (out, lambda.grid = NA, cost = 1, nstep.cv = 80, rmse = NA, 
    offset = 0, trmin = NA, trmax = NA, verbose = FALSE, tol = 1e-05) 
{
    shat.pure.error <- out$shat.pure.error
    pure.ss <- out$pure.ss
    nt <- 2
    np <- out$np
    N <- out$N
    out$cost <- cost
    out$offset <- offset

# find good end points for lambda coarse grid. 

    if (is.na(trmin)) 
        trmin <- 2.05
    if (is.na(trmax)) 
        trmax <- out$np * 0.95

    if (verbose) {
        cat( "trmin and trmax", fill=TRUE)
        cat(trmin, trmax, fill = TRUE)
    }

    if (is.na(lambda.grid[1])) {
        l2 <- sreg.df.to.lambda(trmax, out$xM, out$weightsM)
        l1 <- sreg.df.to.lambda(trmin, out$xM, out$weightsM)
        lambda.grid <- exp(seq(log(l2), log(l1), , nstep.cv))
    }
    if (verbose) {
        cat( "endpoints of coarse lamdba grid", fill=TRUE)
        cat(l1, l2, fill = TRUE)
    }

# build up table of coarse grid serach results for lambda
# in the matrix gcv.grid

    nl <- length(lambda.grid)
    V <- V.model <- V.one <- trA <- MSE <- RSS.model <- rep(NA, 
        nl)

#   loop through lambda's and compute various quantities related to 
#   lambda and the fitted spline. 
    for (k in 1:nl) {
        temp <- sreg.fit(lambda.grid[k], out, verbose=verbose)
        RSS.model[k] <- temp$rss
        trA[k] <- temp$trace
        V[k] <- temp$gcv
        V.one[k] <- temp$gcv.one
        V.model[k] <- temp$gcv.model
    }

# adjustments to columns of gcv.grid
    RSS <- RSS.model + pure.ss
    shat <- sqrt(RSS/(N - trA))
    gcv.grid <- cbind(lambda.grid, trA, V, V.one, V.model, shat)
    dimnames(gcv.grid) <- list(NULL, c("lambda", "trA", "GCV", 
        "GCV.one", "GCV.model", "shat"))

    if (verbose) {
        cat( "Results of coarse grid search", fill=TRUE)
        print(gcv.grid)
    }

    lambda.est <- matrix(ncol = 4, nrow = 5, dimnames = list(c("GCV", 
        "GCV.model", "GCV.one", "RMSE", "pure error"), c("lambda", 
        "trA", "GCV", "shat")))

# now do various refinements for different flavors of finding 
# a good value for lambda the smoothing parameter

##### traditional leave-one-out
    lambda.est[1, 1] <- Krig.find.gcvmin(out, lambda.grid, gcv.grid[, 
        "GCV"], sreg.fgcv, tol = tol, verbose = verbose)
    if (verbose) {
        cat("leave one out GCV lambda", fill=TRUE)
        cat(lambda.est[1, 1], fill = TRUE)
    }

##### using GCV criterion adjusting for replicates
    if (!is.na(shat.pure.error)) {
        lambda.est[2, 1] <- Krig.find.gcvmin(out, lambda.grid, 
            gcv.grid[, "GCV.model"], sreg.fgcv.model, tol = tol, 
            verbose = verbose)
        if (verbose) {
            cat("results of GCV.model search", fill=TRUE) 
            cat(lambda.est[2, 1], fill = TRUE)
        }
    }
    lambda.est[3, 1] <- Krig.find.gcvmin(out, lambda.grid, gcv.grid[, 
        "GCV.one"], sreg.fgcv.one, tol = tol, verbose = verbose)

##### matching an external value of RMSE
    lambda.rmse <- NA
    lambda.pure.error <- NA
    if (!is.na(rmse)) {
        guess <- max(gcv.grid[gcv.grid[, "shat"] < rmse, "lambda"])

        if (verbose) {
            cat("trying to matching with RMSE", fill=TRUE) 
            cat("rmse", rmse, "guess", guess, fill=TRUE)
        }
   
     if (!is.na(guess)) {
            lambda.rmse <- find.upcross(sreg.fs2hat, out, upcross.level = rmse^2, 
                guess = guess, tol = tol * rmse^2)
            lambda.est[4, 1] <- lambda.rmse
   
        }
        else {
            warning("Value of rmse is outside possible range")
        }
    }

##### matching  sigma estimated from the replicates. 
    if (!is.na(shat.pure.error)) {
#       all shats smaller than pure error estimate
        guess <- gcv.grid[, "shat"] < shat.pure.error 

#       set to NA a bad guess!
        if( any( guess) ) { 
             guess<- max(  gcv.grid[guess, "lambda"]) }
        else{
           guess<- NA}

        if (verbose) {
            cat("#### trying to matching with sigma from pure error", 
fill=TRUE)
            cat("shat.pure", shat.pure.error, "guess", guess, fill=TRUE)
        }

        if (!is.na(guess)) {
            lambda.pure.error <- find.upcross(sreg.fs2hat, out, 
                upcross.level = shat.pure.error^2, guess = guess, 
                tol = tol * shat.pure.error^2)
            lambda.est[5, 1] <- lambda.pure.error
            cat("results of matching with pure error sigma", fill=TRUE) 
        }
        else {
            warning("Value of pure error estimate   
                             is outside possible range")
        }
    }
    if (verbose) {
        cat("All forms of estimated lambdas so far", fill=TRUE)
        print(lambda.est)
    }
    for (k in 1:5) {
        lam <- lambda.est[k, 1]
        if (!is.na(lam)) {
            temp <- sreg.fit(lam, out)
            lambda.est[k, 2] <- temp$trace
            if ((k == 1) | (k > 3)) {
                lambda.est[k, 3] <- temp$gcv
            }
            if (k == 2) {
                lambda.est[k, 3] <- temp$gcv.model
            }
            if (k == 3) {
                lambda.est[k, 3] <- temp$gcv.one
            }
            lambda.est[k, 4] <- temp$shat
        }
    }
    list(gcv.grid = gcv.grid, lambda.est = lambda.est)
}

# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"sreg" <-
function(x, y, lam = NA, df = NA, offset = 0, 
        weights = rep(1, length(x)), cost = 1, 
	nstep.cv = 80, find.diagA = TRUE, trmin = 2.01, trmax = 
	length(unique(x)) * 0.95, lammin = NA, lammax = NA,
	verbose = FALSE, do.cv = TRUE, method = "GCV", rmse = NA, lambda = NA,
        na.rm=TRUE)
{
	call <- match.call()
	out <- list()
	out$call <- match.call()
	class(out) <- c("sreg")

	out$cost <- cost
	out$offset <- offset
	out$method <- method
## some obscure components so that some of the Krig functions
## work with out
# size of "null space" 
	out$nt <- 2
        out$knots<-NULL
#
# various checks on x and y including looking for NAs.

        out2<- Krig.check.xY( x,y,NULL, weights, na.rm, verbose=verbose)
        out<- c( out,out2)
        
	## find duplicate rows of the x vector
        ## unique x values are now in out$xM and the means of 
        ## y are in out$yM.
        out<- Krig.replicates( out, verbose=verbose)
        out<- c( out,out2)

# number of unique locations
        out$np<- length( out$yM)

        if( verbose) { print(out)}
	#

        # sorted unique values for prediction to make line plotting quick
        xgrid <- sort(out$xM)
	out$trace <- NA
        
	#
	# figure out if the GCV function should be minimized
	# and that value of lambda should be used for the estimate
	#
	if(!is.na(lambda[1])) {
		lam <- lambda
	}
	if(is.na(lam[1]) & is.na(df[1])) {
		do.cv <- TRUE
	}
	else {
		do.cv <- FALSE
	}

        # find lambda's if df's are given
        #
	if(!is.na(df[1])) {
		lam <- rep(0, length(df))
		for(k in 1:length(df)) {
			lam[k] <- sreg.df.to.lambda(
                                df[k], out$xM, out$weightsM)
		}
	}

############################
        if( verbose) { 
           cat("lambda grid",fill=TRUE)
           print(lam)}
############################

	if(do.cv) {
		a <- gcv.sreg(out, lambda = lam, cost = cost, offset = offset,
			nstep.cv = nstep.cv, verbose = verbose, trmin = trmin,
			trmax = trmax,rmse = rmse)

		# if the spline should be evaluated at the GCV solution 
                # wipe out lam grid
		# and just use GCV lambda.
		out$gcv.grid <- a$gcv.grid
		out$lambda.est <- a$lambda.est
		# 
		# save  GCV estimate if that is what is needed
		lam <- a$lambda.est[method, "lambda"]
		out$shat.GCV <- a$lambda.est[method, "shat"]
	}

# now evaluate spline at lambda either from specified grid or GCV value.

	b <- list()
	# lam can either be  a grid or just the GCV value 
	NL <- length(lam)
	NG <- length(xgrid)
	h <- log(lam)
	residuals <- matrix(NA, ncol = NL, nrow = length(out$y))
	predicted <- matrix(NA, ncol = NL, nrow = NG)
        trace<- rep( NA, NL)
	job <- as.integer(c(0, 3, 0))

	if(find.diagA) {
		diagA <- matrix(NA, ncol = NL, nrow = out$np)
		# add switch to find diag of A. 
		job <- as.integer(c(3, 3, 0))
	}

	for(k in 1:NL) {
		#
		# call cubic spline FORTRAN, this is nasty looking but fast.
                # note lambda is passed in log scale.
                # what the routine does is controlled by array job 
		# spline solution evaluated at xgrid
		# 
		b <- .Fortran("css",
			h = as.double(h[k]),
			npoint = as.integer(out$np),
			x = as.double(out$xM),
			y = as.double(out$yM),
			wt = as.double(1/sqrt(out$weightsM)),
			sy = as.double(rep(0, out$np)),
			trace = as.double(0),
			diag = as.double(c(cost, offset, rep(0, (out$np - 2)))),
       			cv = as.double(0),
			ngrid = as.integer(NG),
			xg = as.double(xgrid),
			yg = as.double(rep(0, NG)),
			job = as.integer(job),
			ideriv = as.integer(0),
			ierr = as.integer(0), PACKAGE="fields")
		if(find.diagA) {
		diagA[, k] <- b$diag}

# note distinction between yM and y, xM and x
# these are residuals at all data point locations not just the
# unique set. 
	        trace[k]<- b$trace
        	residuals[, k] <- out$y - splint(out$xM, b$sy, out$x)
		predicted[, k] <- b$yg }
	
        out$call <- call
	out$lambda <- lam
	out$do.cv <- do.cv

	out$residuals <- residuals
        out$trace<- trace
	out$fitted.values <- out$y - residuals
	out$predicted <- list(x = xgrid, y = predicted)

	if(length(lambda[1]) == 1) {
		out$eff.df <- out$trace[1]
	}

	if(find.diagA) {
		out$diagA <- diagA
	}

        class(out)<- "sreg"

	return(out)
}

# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"sreg.df.to.lambda" <-
function(df, x, wt, guess = 1, tol = 1.0000000000000001e-05)
{

	if(is.na(df))
		return(NA)
	n <- length(unique(x))
	info <- list(x = x, wt = wt, df = df)
	if(df > n) {
		warning(" df too large to match a lambda value")
		return(NA)
	}
	h1 <- log(guess)

	########## find upper lambda
	for(k in 1:25) {
		tr <- sreg.trace(h1, info)
		if(tr <= df)
			break
		h1 <- h1 + 1.5
	}

	########## find lower lambda
	h2 <- log(guess)
	for(k in 1:25) {
		tr <- sreg.trace(h2, info)
		if(tr >= df)
			break
		h2 <- h2 - 1.5
	}
	out <- bisection.search(h1, h2, sreg.fdf, tol = tol, f.extra = info)$
		x
	exp(out)
}
# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"sreg.fdf" <-
function(h, info)
{
	sreg.trace(h, info) - info$df
}
# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"sreg.fgcv" <-
function(lam, obj)
{
	sreg.fit(lam, obj)$gcv
}
# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"sreg.fgcv.model" <-
function(lam, obj)
{
	sreg.fit(lam, obj)$gcv.model
}
# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"sreg.fgcv.one" <-
function(lam, obj)
{
	sreg.fit(lam, obj)$gcv.one
}
# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"sreg.fit" <-
function(lam, obj, verbose=FALSE)
{
	np <- obj$np
	N <- obj$N
	nt <- 2
	if(is.null(obj$cost)) {
		cost <- 1
	}
	else {
		cost <- obj$cost
	}
	if(is.null(obj$offset)) {
		offset <- 0
	}
	else {
		offset <- obj$offset
	}
	if(is.null(obj$shat.pure.error)) {
		shat.pure.error <- 0
	}
	else {
		shat.pure.error <- obj$shat.pure.error
	}
	if(is.null(obj$pure.ss)) {
		pure.ss <- 0
	}
	else {
		pure.ss <- obj$pure.ss
	}
	#print(np)
	#	NOTE h <- log(lam)
	temp <- .Fortran("css",
		h = as.double(log(lam)),
		npoint = as.integer(np),
		x = as.double(obj$xM),
		y = as.double(obj$yM),
		wt = as.double(sqrt(1/obj$weightsM)),
		sy = as.double(rep(0, np)),
		trace = as.double(0),
		diag = as.double(rep(0, np)),
		cv = as.double(0),
		ngrid = as.integer(0),
		xg = as.double(0),
		yg = as.double(0),
		job = as.integer(c(3, 0, 0)),
		ideriv = as.integer(0),
		ierr = as.integer(0),
                PACKAGE="fields")
	rss <- sum((temp$sy - obj$yM)^2 * obj$weightsM)
	MSE <- rss/np
	if((N - np) > 0) {
		MSE <- MSE + pure.ss/(N - np)
	}
	trA <- temp$trace
	den <- (1 - (cost * (trA - nt - offset) + nt)/np)
	den1 <- (1 - (cost * (trA - nt - offset) + nt)/N)

	# If the denominator is negative then flag this as a bogus case
	# by making the GCV function "infinity"
	#
	shat <- sqrt((rss + pure.ss)/(N - trA))
	GCV <- ifelse(den > 0, MSE/den^2, NA)
	gcv.model <- ifelse(
              (den > 0)&( (N-np)>0), 
              pure.ss/(N - np) + (rss/np)/(den^2), NA)
	gcv.one <- ifelse(den > 0, ((pure.ss + rss)/N)/(den1^2), NA)
	list(trace = trA, gcv = GCV, rss = rss, shat = shat, gcv.model = 
		gcv.model, gcv.one = gcv.one)
}

# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"sreg.fs2hat" <-
function(lam, obj)
{
	sreg.fit(lam, obj)$shat^2
}
# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"sreg.trace" <-
function(h, info)
{

	N <- length(info$x)
	#	h <- log(lam)
	temp <- (.Fortran("css",
		h = as.double(h),
		npoint = as.integer(N),
		x = as.double(info$x),
		y = as.double(rep(0, N)),
		wt = as.double(1/sqrt(info$wt)),
		sy = as.double(rep(0, N)),
		trace = as.double(0),
		diag = as.double(rep(0, N)),
		cv = as.double(0),
		ngrid = as.integer(0),
		xg = as.double(0),
		yg = as.double(0),
		job = as.integer(c(3, 0, 0)),
		ideriv = as.integer(0),
		ierr = as.integer(0),PACKAGE="fields")$trace)
	return(temp)
}
