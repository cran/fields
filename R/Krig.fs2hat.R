# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"Krig.fs2hat" <-
function(lam, obj)
{
	lD <- obj$matrices$D * lam
	RSS <- obj$pure.ss + sum(((obj$matrices$u * lD)/(1 + lD))^2)
	#	print(RSS)
	#	trA <- sum(1/(1 + lD)) + obj$offset
	den <- obj$N - (sum(1/(1 + lD)) + obj$offset)
	if(den < 0) {
		return(NA)
	}
	else {
		RSS/(den)
	}
}
