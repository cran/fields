"dup" <-
function (dat) 
{
    dat <- match(dat, unique(dat))
    id <- order(dat)
    look <- c(FALSE, ifelse(diff(sort(dat)) == 0, TRUE, FALSE))
    dat[id] <- look
    as.logical(dat)
}
