"exp.cov.simple" <-
function (x1, x2, theta = 1) 
{
    exp(-rdist(x1, x2)/theta)
}
