# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
mKrig.MLE <- function(x, y, weights = rep(1, nrow(x)), cov.fun="stationary.cov", cov.args = NULL, 
                      Z = NULL, par.grid = NULL, lambda = NULL, lambda.profile = TRUE, 
                      verbose = FALSE, relative.tolerance = 1e-04, ...) {
  
  #check which optimization options the covariance function supports
  supportsDistMat = supportsArg(cov.fun, "distMat")
  
  #precompute distance matrix if possible so it only needs to be computed once
  if(supportsDistMat) {
    
    #Get distance function and arguments if available.  Otherwise use 'dist' function
    #to compute upper triangle of distance matrix
    #
    Dist.fun= c(cov.args, list(...))$Distance
    Dist.args=c(cov.args, list(...))$Dist.args
    
    if(is.null(Dist.fun))
      Dist.fun = "dist"
    
    distMat = do.call(Dist.fun, c(list(x), Dist.args))
  }
  
  # mKrig.args has all the arguments needed to call mKrig except lambda and cov.args
  if(supportsDistMat)
    cov.args = c(cov.args, list(distMat=distMat, onlyUpper=TRUE))
  mKrig.args <- c(list(x = x, y = y, weights = weights, Z = Z, cov.fun=cov.fun), 
                  list(...))
  
  lnProfileLike.max <- -1e+20
  
  # find NG --  number of parameters to try
  par.grid <- data.frame(par.grid)
  if (nrow(par.grid) == 0) {
    if (is.null(lambda)) {
      NG <- 1
    }
    else {
      NG <- length(lambda)
    }
  }
  else {
    NG <- nrow(par.grid)
  }
  
  # output matrix to summarize results
  summary <- matrix(NA, nrow = NG, ncol = 8)
  dimnames(summary) <- list(NULL, c("EffDf", "lnProfLike", 
                                    "GCV", "sigma.MLE", "rho.MLE", "llambda.MLE", "counts eval", 
                                    "counts grad"))
  lambda.best <- NA
  
  # default for lambda is 1.0 for first value and exp(llambda.opt) for subsequent ones
  # this is controlled by NAs for lambda starting values.
  if (is.null(lambda)) {
    lambda <- rep(NA, NG)
  }
  
  # default starting value for lambda is 1 or log lambda is 0
  llambda.opt <- 0
  optim.counts <- c(NA, NA)
  lnLike.eval <- list()
  
  # Define the objective function as a tricksy call to mKrig
  # if Y is a matrix of replicated data sets use the log likelihood for the complete data sets
  temp.fn <- function(x) {
    # NOTE: FULL refers to estimates collapsed across the replicates if Y is a matrix
    # assign to hold only a few components returned by mKrig
    hold <- do.call("mKrig", c(mKrig.args, list(find.trA = FALSE, lambda = exp(x), 
                                                cov.args=c(cov.args.temp, cov.args)))
                    )[c("lambda.fixed", "rho.MLE.FULL", "sigma.MLE.FULL", "lnProfileLike.FULL")]
    
    # add this evalution to an  object (i.e. here a matrix) in the calling frame
    temp.eval <- get("capture.evaluations")
    assign("capture.evaluations", rbind(temp.eval, unlist(hold)), 
           envir = capture.env)
    return(hold$lnProfileLike.FULL)
  }
  #
  # begin loop over covariance arguments
  for (k in 1:NG) {
    llambda.start <- ifelse(is.na(lambda[k]), llambda.opt, log(lambda[k]))
    
    # list of covariance arguments from par.grid with right names (some R arcania!)
    # note that this only works because 1) temp.fn will search in this frame for this object
    # par.grid has been coerced to a data frame so one has a concept of a row subscript.
    cov.args.temp <- as.list(par.grid[k, ])
    names(cov.args.temp) <- names(par.grid)
    
    #optimize over lambda if lambda.profile is TRUE
    if (lambda.profile) {
      # set up matrix to store evaluations from within optim
      capture.evaluations <- matrix(NA, ncol = 4, nrow = 1, 
                                    dimnames = list(NULL, c("lambda", "rho.MLE", 
                                                            "sigma.MLE", "lnProfileLike.FULL")))
      capture.env <- environment()
      
      # call to optim
      look <- optim(llambda.start, temp.fn, method = "BFGS", 
                    control = list(fnscale = -1, parscale = 0.1, 
                                   ndeps = 0.05, reltol = relative.tolerance))
      llambda.opt <- look$par
      optim.counts <- look$counts
      
      # call to 1-d search
      #            opt.summary <- optimize(temp.fn, interval= llambda.start + c(-8,8), maximum=TRUE)
      #            llambda.opt <- opt.summary$maximum
      #            optim.counts<- c(nrow(capture.evaluations)-1, NA)
      
      # accumulate the new matrix of lnlambda and ln likelihoods (omitting first row of NAs)
      lnLike.eval <- c(lnLike.eval, list(capture.evaluations[-1, ]))
    }
    else {
      # no refinement for lambda so just save the the 'start' value as final one.
      llambda.opt <- llambda.start
    }
    
    # final fit at optimal value (or starting value if not refinement/maximization for lambda)
    obj <- do.call("mKrig", c(mKrig.args, list(lambda = exp(llambda.opt), 
                                               cov.args=c(cov.args.temp, cov.args))))
    if (obj$lnProfileLike.FULL > lnProfileLike.max) {
      lnProfileLike.max <- obj$lnProfileLike.FULL
      cov.args.MLE <- cov.args.temp
      lambda.best <- exp(llambda.opt)
    }
    
    # save results of the kth covariance model evaluation
    summary[k, 1:8] <- c(obj$eff.df, obj$lnProfileLike.FULL, 
                         obj$GCV, obj$sigma.MLE.FULL, obj$rho.MLE.FULL, llambda.opt, 
                         optim.counts)
    if (verbose) {
      cat("Summary: ", k, summary[k, 1:8], fill = TRUE)
    }
  }
  return(list(summary = summary, par.grid = par.grid, cov.args.MLE = cov.args.MLE, 
              mKrig.args = list(...), lambda.best = lambda.best, lambda.MLE = lambda.best, 
              call = match.call(), lnLike.eval = lnLike.eval))
}
