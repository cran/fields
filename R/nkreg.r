"nkreg" <-
function(data.x, data.y, bandwidth, n.points=50, grid=NULL,
grid.list=NULL)
{
	# who wants to keep typing bandwidth!
	h <- bandwidth
	# reformat 1-d vector as column matrix
	if(data.class(data.x) == "data.frame") data.x <- as.matrix(data.x)
	if(!is.matrix(data.x))
		m <- 1
	else m <- ncol(data.x)
#
# special case for 2-d 
#
if( is.null( grid.list)& is.null( grid) &  (m==2)){
xr<- range( data.x[,1])
yr<- range( data.x[,2])
grid.list<- list( 
seq( xr[1], xr[2],,n.points),
seq( yr[1], yr[2],,n.points))
class( grid.list)<- "grid.list"
}
if( !is.null(grid.list)){
grid<- make.surface.grid( grid.list)
}
# special case for 1-d
	if(is.null(grid)) {
		if((m == 1))
				grid <- seq(min(data.x), max(data.x),  , 
					n.points)
		else grid <- data.x
	}
	if(data.class(grid) == "data.frame")
		grid <- as.matrix(grid)
	if(m > 1) {
		if(m != ncol(grid))
			stop("dimension of data.x and grid do not agree")
		p <- nrow(grid)
		n <- nrow(data.x)
	} else {
		p <- length(grid)
		n <- length(data.y)
	}
	nh <- length(h)
	f <- matrix(rep(-99, p * nh), ncol = nh)
	out <- list()
	for(k in 1:nh) {
		f[,k] <- .Fortran("nkreg",
			as.double(h[k]),
			as.integer(n),
			as.integer(m),
			as.double(data.x),
			as.double(data.y),
			as.integer(p),
			as.double(grid),
			estimate=as.double(rep(-99, p)))$estimate
	}
	if(nh == 1)
		f <- c(f)
return(
	list(x = grid, y = f, h = h, grid.list=grid.list))
}
