"nkden" <-
function(data, bandwidth, n.points, grid)
{
	if(!is.loaded(symbol.For("nkden"))) {
		temp <- dyn.load(paste(FUNFITS.BIN, "funfits.o", sep = ""),
			2)
	}
	# who wants to keep typing bandwidth!]
	if(data.class(data) == "data.frame") data <- as.matrix(data)
	h <- bandwidth
	# reformat 1-d vector as column matrix
	if(!is.matrix(data)) m <- 1 else m <- ncol(data)
	if(missing(grid)) {
		if((m == 1))
			if(!missing(n.points))
				grid <- seq(min(data), max(data),  , n.points)
			else grid <- sort(data)
		else {
			grid <- data
		}
	}
	if(data.class(grid) == "data.frame")
		grid <- as.matrix(grid)
	if(m > 1) {
		if(m != ncol(grid))
			stop("dimension of data and grid do not agree")
		p <- nrow(grid)
		n <- nrow(data)
	} else {
		p <- length(grid)
		n <- length(data)
	}
	nh <- length(h)
	f <- matrix(rep(-99, p * nh), ncol = nh)
	out <- list()
	for(k in 1:nh) {
		out <- .Fortran("nkden",
			as.double(h[k]),
			as.integer(n),
			as.integer(m),
			as.double(data),
			as.integer(p),
			as.double(grid),
			as.double(rep(-99, p)))
		f[, k] <- out[[7]]
	}
	if(nh == 1)
		f <- c(f)
	list(x = grid, y = f, h = h)
}
