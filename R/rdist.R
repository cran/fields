# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"rdist" <- function(x1, x2 = NULL, compact = FALSE) {
	if (!is.matrix(x1)) {
		x1 <- as.matrix(x1)
	}
  if (is.null(x2)) {
		storage.mode(x1) <- "double"
		if (compact)
			  return(dist(x1))
		else
			return(.Call("RdistC", x1, x1, PACKAGE = "fields"))
	} else {
	  if (!is.matrix(x2)) {
			x2 <- as.matrix(x2)
		}
		storage.mode(x1) <- "double"
		storage.mode(x2) <- "double"
		return(.Call("RdistC", x1, x2, PACKAGE = "fields"))
	}

}
