"print.summary.Krig" <-
function(x, ...)
{
	digits <- x$digits
	c1 <- "Number of Observations:"
	c2 <- x$num.observation
	c1 <- c(c1, "Number of unique points:")
	c2 <- c(c2, x$num.uniq)
	c1 <- c(c1, "Degree of polynomial null space ( base model):")
	c2 <- c(c2, x$m - 1)
	c1 <- c(c1, "Number of parameters in the null space")
	c2 <- c(c2, x$nt)
	c1 <- c(c1, "Effective degrees of freedom:")
	c2 <- c(c2, format(round(x$enp, 1)))
	c1 <- c(c1, "Residual degrees of freedom:")
	c2 <- c(c2, format(round(x$num.observation - x$enp, 1)))
	c1 <- c(c1, "MLE sigma ")
	c2 <- c(c2, format(signif(x$shat.MLE, digits)))
	c1 <- c(c1, "GCV est. sigma ")
	c2 <- c(c2, format(signif(x$shat.GCV, digits)))
	if(!is.na(x$shat.pure.error)) {
		c1 <- c(c1, "Pure error sigma")
		c2 <- c(c2, format(signif(x$shat.pure.error, digits)))
	}
	c1 <- c(c1, "MLE rho ")
	c2 <- c(c2, format(signif(x$rhohat, digits)))
	c1 <- c(c1, "Scale used for covariance (rho)")
	c2 <- c(c2, signif(x$rho, digits))
	c1 <- c(c1, "Scale used for nugget (sigma^2)")
	c2 <- c(c2, signif(x$sigma2, digits))
	c1 <- c(c1, "lambda (sigma2/rho)")
	c2 <- c(c2, signif(x$lambda, digits))
	#	c1 <- c(c1, "Cost in GCV")
	#	c2 <- c(c2, format(round(x$cost, 2)))
	#	c1 <- c(c1, "GCV Minimum")
	#	c2 <- c(c2, format(signif(x$gcvmin, digits)))
	sum <- cbind(c1, c2)
	dimnames(sum) <- list(rep("", dim(sum)[1]), rep("", dim(sum)[2]))
	res.quantile <- x$res.quantile
	names(res.quantile) <- c("min", "1st Q", "median", "3rd Q", "max")
	###
	###
	###
	cat("CALL:\n")
	dput(x$call)
	print(sum, quote = FALSE)
	cat("\n")
	cat("RESIDUAL SUMMARY:", fill = TRUE)
	print(signif(res.quantile, digits))
	cat("\n")
	cat("COVARIANCE MODEL:", x$cov.function, fill = TRUE)
	if((x$correlation.model)) {
		cat(" A correlation model was fit: Y is standardized\nbefore\nspatial estimate is found",
			fill = TRUE)
	}
	if(x$knot.model) {
		cat(" Knot model: ", x$np, 
			" knots supplied to define basis functions", fill = TRUE)
	}
	cat("\n")
	cat("DETAILS ON SMOOTHING PARAMETER:", fill = TRUE)
	cat(" Method used:  ", x$method, "   Cost: ", x$cost, fill = TRUE)
	#	cat(" Stats on this choice of lambda", fill = TRUE)
	print(x$sum.gcv.lambda, digits = digits)
	cat("\n")
	cat(" Summary of estimates for lambda", fill = TRUE)
	print(x$lambda.est, digits = x$digits)
	invisible(x)
}
