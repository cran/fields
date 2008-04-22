# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"Krig.fdf" <-
function(llam, info)
{
	sum(1/(1 + exp(llam) * info$D)) - info$df
}
