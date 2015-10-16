spatialProcess <- function(x, y, cov.function = "stationary.cov", 
	cov.args = list(Covariance = "Matern", smoothness = 1),
	 ngrid=10, theta.grid = NULL, ...) {
	MLEfit <- MLESpatialProcess(x, y, cov.function = cov.function, 
		cov.args = cov.args, ngrid=ngrid, theta.grid = theta.grid,
		 ...)
# now fit spatial model with MLE for theta (range parameter)
# reestimate the other parameters for simplicity to get the complete Krig
# object.		
	obj <- Krig(x, y, cov.function = cov.function, cov.args = cov.args, 
		theta = MLEfit$pars[1],  
		method = "REML", give.warnings=TRUE, 
		...)
	obj <- c(obj, MLEfit)
	obj$theta.MLE<- MLEfit$pars[1]
# replace call with this top level one
    obj$call<- match.call()	
	class(obj) <- c( "spatialProcess","Krig")
 
	return(obj)
}
