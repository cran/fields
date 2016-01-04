# fields, Tools for spatial data
# Copyright 2015, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
test.for.zero <- function(xtest, xtrue, tol = 1e-08, 
    relative = TRUE, tag = NULL) {
    denom <- ifelse(relative, mean(abs(c(xtrue))), 1)
    test.value <- sum(abs(c(xtest) - c(xtrue)))/denom
    if (!is.null(tag)) {
      cat("Testing: ", tag, fill = TRUE)
    }
    if (test.value < tol) {
            cat("PASSED test at tolerance ", tol, fill = TRUE)
    }
    else {

          cat("FAILED test value = ", test.value, " at tolerance ", 
            tol)
# generate an "error" to signal failed test      
      if (exists("test.for.zero.flag")) {
        stop()
      }
    }
}
