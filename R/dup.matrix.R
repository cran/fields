"dup.matrix" <-
function (mat) 
{
    mat <- as.matrix(mat)
    nc <- ncol(mat)
    temp <- matrix(match(c(mat), unique(c(mat))), ncol = nc)
    temp2 <- format(temp[, 1])
    if (nc > 1) {
        for (k in 2:nc) {
            temp2 <- paste(temp2, temp[, k], sep = "X")
        }
    }
    dup(temp2)
}
