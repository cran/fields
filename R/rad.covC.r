"rad.covC" <-
function (x1, x2, C, p = 1, with.log = TRUE, with.constant = TRUE) 
{
    if (!is.loaded(symbol.For("radbas"))) {
        temp <- dyn.load(paste(FIELDS.BIN, "fields.o", sep = ""), 
            2)
    }
    if (!is.matrix(x1)) {
        x1 <- as.matrix(x1)
    }
    if (!is.matrix(x2)) {
        x2 <- as.matrix(x2)
    }
    d <- ncol(x1)
    n1 <- nrow(x1)
    n2 <- nrow(x2)
    par <- c(p/2, ifelse(d%%2 == 0, 1, 0))
    if (!with.log) 
        par[2] <- 0
    temp <- .Fortran("multrb", nd = as.integer(d), x1 = as.double(x1), 
        n1 = as.integer(n1), x2 = as.double(x2), n2 = as.integer(n2), 
        par = as.double(par), c = as.double(C), h = as.double(rep(0, 
            n1)), work = as.double(rep(0, n2)))
    if (with.constant) {
        m <- (d + p)/2
        Amd <- radbas.constant(m, d)
    }
    else {
        Amd <- 1
    }
    Amd * temp$h
}
