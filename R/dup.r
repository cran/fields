"dup" <-
function (dat) 
{
    dat <- match(dat, unique(dat))
    id <- order(dat)
    look <- c(F, ifelse(diff(sort(dat)) == 0, T, F))
    dat[id] <- look
    as.logical(dat)
}
