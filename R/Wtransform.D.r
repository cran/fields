"Wtransform.D" <-
function (nx, ny, weights = c(1), cut.min = 8, details = F) 
{
    NN <- ny
    MM <- nx
    while ((NN > cut.min) & (MM > cut.min)) {
        NN <- NN/2
        MM <- MM/2
    }
    temp <- matrix(0, nrow = nx, ncol = ny)
    nn <- NN
    mm <- MM
    temp[1:mm, 1:nn] <- weights[1]
    cat(mm, nn, fill = T)
    jj <- 2
    nn <- nn * 2
    mm <- mm * 2
    while ((nn <= ny) & (mm <= nx) & (jj <= length(weights))) {
        cat(mm, nn, fill = T)
        if (!details) {
            temp[(mm/2 + 1):mm, 1:(nn/2)] <- weights[jj]
            temp[(mm/2 + 1):mm, (nn/2 + 1):nn] <- weights[jj]
            temp[1:(mm/2), (nn/2 + 1):nn] <- weights[jj]
        }
        else {
            temp[(mm/2 + 1):mm, 1:(nn/2)] <- weights[jj]
            jj <- jj + 1
            temp[(mm/2 + 1):mm, (nn/2 + 1):nn] <- weights[jj]
            jj <- jj + 1
            temp[1:(mm/2), (nn/2 + 1):nn] <- weights[jj]
        }
        jj <- jj + 1
        nn <- nn * 2
        mm <- mm * 2
    }
    max.m <- mm/2
    max.n <- nn/2
    return(D = temp, max.n, max.m)
}
