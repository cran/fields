# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
double.exp <- function(x) {
    # double exponential weight function
    0.5 * exp(-abs(x))
}
