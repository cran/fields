"WQS.basis" <-
function (N, cut.n = 8) 
{
    x <- diag(1, N)
    nrow <- nrow(x)
    ncol <- ncol(x)
    nn <- nrow
    temp <- x
    nn <- cut.n * 2
    while (nn <= nrow) {
        temp[1:nn, ] <- WQSi(temp[1:nn, ])
        nn <- nn * 2
    }
    return(temp)
}
