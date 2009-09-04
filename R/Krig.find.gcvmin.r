# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"Krig.find.gcvmin" <- function(info, lambda.grid, 
    gcv, gcv.fun, tol, verbose = FALSE, give.warnings = TRUE) {
    ind <- !is.na(gcv)
    lambda.grid <- lambda.grid[ind]
    gcv <- gcv[ind]
    nstep.cv <- length(lambda.grid)
    il <- order(gcv)[1]
    lambda.gcv <- lambda.grid[il]
    gcv.raw <- min(gcv)
    if (verbose) {
        cat("#### Call for refined search using", as.character(substitute(gcv.fun)), 
            fill = TRUE)
        cat("Results of coarse search lambda and GCV:", lambda.grid[il], 
            gcv.raw, fill = TRUE)
    }
    if ((il > 1) & (il < nstep.cv)) {
        out <- golden.section.search(lambda.grid[il - 1], lambda.grid[il], 
            lambda.grid[il + 1], gcv.fun, f.extra = info, tol = tol * 
                gcv.raw)
        return(out$x)
    }
    else {
        if (give.warnings) {
            warning("GCV search gives a minimum at the endpoints of the\ngrid search")
        }
        return(lambda.gcv)
    }
}
