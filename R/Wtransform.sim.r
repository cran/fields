"Wtransform.sim" <-
function (D, cut.min = 16) 
{
    NN <- dim(D)
    temp <- matrix(rnorm(NN[1] * NN[2]), NN[1], NN[2]) * sqrt(D)
    Wtransform.image(temp, inv = T, cut.min = cut.min)
}
