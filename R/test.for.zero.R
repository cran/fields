# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
test.for.zero <- function(xtest, xtrue, tol = 1e-08, 
    relative = TRUE, tag = NULL) {
    denom <- ifelse(relative, mean(abs(c(xtrue))), 1)
    test.value <- sum(abs(c(xtest) - c(xtrue)))/denom
    if (test.value < tol) {
        if (exists("test.for.zero.flag")) {
            if (!is.null(tag)) {
                cat("testing: ", tag, fill = TRUE)
            }
            cat("    PASSED test at tolerance ", tol, fill = TRUE)
        }
    }
    else {
        if (!is.null(tag)) {
            cat("testing: ", tag, fill = TRUE)
        }
        cat("FAILED test value = ", test.value, " at tolerance ", 
            tol, fill = TRUE)
    }
}
