"in.poly" <-
function(xd, xp, convex.hull = FALSE)
{
	if(!is.loaded(symbol.For("inpoly"))) {
		temp <- dyn.load(paste(FIELDS.BIN, "fields.o", sep = ""), 2)
	}
	if(convex.hull) {
		xp <- xp[chull(xp),  ]
	}
	nd <- nrow(xd)
	np <- as.integer(nrow(xp))
	.Fortran("inpoly",
		nd = as.integer(nd),
		as.single(xd[, 1]),
		as.single(xd[, 2]),
		np = np,
		as.single(xp[, 1]),
		as.single(xp[, 2]),
		ind = as.integer(rep(-1, nd)))$ind
}
