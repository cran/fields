"predict.sreg" <-
function(out, x, derivative = 0, model = 1)
{
	if(missing(x))
		x <- out$xraw
	c(splint(out$predicted$x, out$predicted$y[, model], x, derivative = 
		derivative))
}
