# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"cat.to.list" <-
function (x, a) 
{
    a <- as.character(a)
    label <- unique(a)
    out <- as.list(1:length(label))
    names(out) <- label
    for (k in 1:length(label)) {
        out[[k]] <- x[label[k] == a]
        if (length(out[[k]]) == 0) 
            out[[k]] <- NA
    }
    out
}
