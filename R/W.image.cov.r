"W.image.cov" <-
function (ind1, ind2, Y, cov.obj = NULL, setup = F, m = 64, n = 64, 
    ...) 
{
    if (is.null(cov.obj)) {
        cov.obj <- setup.W.image.cov(m, n, ...)
    }
    if (setup) {
        return(cov.obj)
    }
    cut.min <- cov.obj$cut.min
    if (missing(ind1)) {
        Wtransform.image(cov.obj$D * Wtransform.image(Y, inv = T, 
            transpose = T, cut.min = cut.min), inv = T, cut.min = cut.min)
    }
    else {
        temp <- matrix(0, nrow = cov.obj$m, ncol = cov.obj$n)
        if (missing(ind2)) {
            ind2 <- ind1
        }
        temp[ind2] <- Y
        Wtransform.image(cov.obj$D * Wtransform.image(temp, inv = T, 
            transpose = T, cut.min = cut.min), inv = T, cut.min = cut.min)[ind1]
    }
}
