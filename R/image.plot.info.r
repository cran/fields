"image.plot.info" <-
function(...)
{
	temp <- list(...)
	#
	xlim <- NA
	ylim <- NA
	zlim <- NA
	#
	# go through various cases of what these can be
	#
	if(is.list(temp[[1]])) {
		xlim <- range(temp[[1]]$x, na.rm = T)
		ylim <- range(temp[[1]]$y, na.rm = T)
		zlim <- range(temp[[1]]$z, na.rm = T)
	}
	if(is.matrix(temp[[1]])) {
		xlim <- c(0,1)
		ylim <- c(0,1)
		zlim <- range(temp[[1]], na.rm = T)
	}
   if( length( temp)>=3){
	if(is.matrix(temp[[3]])) {
		xlim <- range(temp[[1]], na.rm = T)
		ylim <- range(temp[[2]], na.rm = T)
		zlim <- range(temp[[3]], na.rm = T)
	}
}
	xthere <- match("x", names(temp))
	ythere <- match("y", names(temp))
	zthere <- match("z", names(temp))
	if(!is.na(zthere))
		zlim <- range(temp$z, na.rm = T)
	if(!is.na(xthere))
		xlim <- range(temp$x, na.rm = T)
	if(!is.na(ythere))
		ylim <- range(temp$y, na.rm = T)
	if(!is.null(temp$zlim))
		zlim <- temp$zlim
	if(!is.null(temp$xlim))
		xlim <- temp$xlim
	if(!is.null(temp$ylim))
		ylim <- temp$ylim
	list(xlim = xlim, ylim = ylim, zlim = zlim)
}
