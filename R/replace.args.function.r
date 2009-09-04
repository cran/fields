# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"replace.args.function" <- function(fun, ...) {
    temp <- list(...)
    ntemp <- names(temp)
    fnames <- names(fun)
    if (length(temp) > 0) {
        for (k in 1:length(ntemp)) {
            if (!is.na(match(ntemp[k], fnames))) {
                fun[ntemp[k]] <- temp[ntemp[k]]
            }
        }
    }
    as.function(fun)
}
