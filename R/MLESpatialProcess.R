


MLESpatialProcess <- function(x, y, theta.grid = NULL, cov.function = "stationary.cov", 
	cov.args = list(Covariance = "Matern", smoothness = 1), ngrid = 10, 
	niter = 15, tol = 0.01, Distance = "rdist", 
	nstep.cv= 50, verbose=FALSE,
	...) {
# save arguments for Krig as a list
	
	KrigCallingList <- c(
	# list to pass to the objective function
    info <- list( x = x, 
                  Y = y, 
       cov.function = cov.function,
           cov.args = cov.args, 
		     method = "REML",
		   nstep.cv = nstep.cv,
      give.warnings = FALSE, 
     	   Distance = Distance),
     	    list( ...)
    	   )
# objective function used for grid search and optimization	
	objective.fn <- function(ltheta, returnAll = FALSE) {	
		thetaList<- list( theta=exp( ltheta))
		hold <- do.call(Krig, 
		     c( KrigCallingList,thetaList) )[c("lambda.est", "gcv.grid")]
		minus.lPLike <- hold$lambda.est["REML", 5]
		# add this evalution to an  object
		# (i.e. here a matrix) in the calling frame
        withAddedRow <- rbind(
             get("capture.evaluations", envir = capture.env),
             c(exp(ltheta),  hold$lambda.est["REML",  ])
             )
             if( verbose){
             	print( c(exp(ltheta),  hold$lambda.est["REML",  ]))
             }
		assign("capture.evaluations", withAddedRow, envir = capture.env)
		# return all serach information or just profile like value.        
		if (returnAll) {
			return(hold)
		} else {
			return(minus.lPLike)
		}
	}	
	# set up matrix to cature all  evaluations from within optimization
	capture.evaluations <- matrix(NA, ncol = 7, nrow = 1, dimnames = list(NULL, 
		c("theta", "lambda.MLE", "trA", "GCV", "sigmaGCV", "lnProfileLike","warning")))
	capture.env <- environment()
	#
	# if grid for ranges is missing use  quantiles of pairwise
#distances among data.
#
if (is.null(theta.grid)) {
		# Distances between locations
		pairwiseD<- get(Distance)(x,x)
		pairwiseD<- pairwiseD[col(pairwiseD) > row( pairwiseD) ]
		theta.range <- quantile(pairwiseD , c(0.03, 0.97))
		theta.grid <- seq(theta.range[1], theta.range[2], , ngrid)
	}
	if (length(theta.grid) == 2) {
		theta.grid <- seq(theta.grid[1], theta.grid[2], , ngrid)
	}
	ngrid <- length(theta.grid)
	minus.REML <- rep(NA, ngrid)
# object to hold log likelihood surface as a function of log lambda and 	theta
	logLikelihoodSurface<- list( x = matrix( theta.grid, nrow= ngrid, ncol=nstep.cv ),
	                       y = matrix(NA,nrow=ngrid, ncol= nstep.cv ),
	                       z = matrix( NA, nrow= ngrid,  ncol=nstep.cv))	                       
	# grid search over theta likelihood maximized over
	# lambda (sill and nugget) for each theta
for (j in 1:ngrid) {
	    out<- objective.fn(log(theta.grid[j]),  returnAll=TRUE)
		minus.REML[j] <- out$lambda.est["REML", 5]
		logLikelihoodSurface$y[j,]<- log( out$gcv.grid[,1])
		logLikelihoodSurface$z[j,]<- -1*out$gcv.grid[,7]
	}
#	temp <- cbind(theta.grid, -minus.REML)
#	dimnames(temp) <- list(NULL, c("theta", "logProfileLike"))
# best point for theta from grid search
# NOTE: due to past convenion in Krig  - log likelihood is computed
# so this quantity is _minimized_
	IMIN <- which.min( minus.REML)
	if (IMIN == 1 | IMIN == ngrid) {
		warning("theta (range parameter) at end of search interval:", fill = TRUE)
		theta.MLE <- theta.grid[IMIN]
	    REML.MLE<- -1* minus.REML[IMIN]
	} else {
		# starting interval  for  1-d  optimization
		lstart <- log(theta.grid)[c(IMIN-1, IMIN+1) ]
		out <- optimize(f = objective.fn, interval = c(lstart[1], lstart[2]), maximum = FALSE)
		theta.MLE <- exp(out$minimum)
		REML.MLE <- -1 * out$objective
	}
	eval.grid <- capture.evaluations[-1, ]
	eval.grid[, 6] <- -1 * eval.grid[, 6]
# sort on theta!
	ind<- order( eval.grid[,1])
	eval.grid<- eval.grid[ind,]
#		
    thetaList<- list( theta=theta.MLE)
		hold <- do.call(Krig, c( KrigCallingList, thetaList) )[c("lambda.est", "rho.MLE", "shat.MLE")]
		pars<- c( theta.MLE, hold$lambda.est[6,1], hold$rho.MLE, hold$shat.MLE)
		names( pars) <- c( "theta", "lambda", "rho", "sigma")
		return(list(
		          pars = pars,
		 logLikelihood = REML.MLE,
		     eval.grid = eval.grid,
  logLikelihoodSurface = logLikelihoodSurface,
            lambda.est = hold$lambda.est,
		          call = match.call() )
		       )   		     
}
