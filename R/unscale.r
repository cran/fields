"unscale" <-
function (x, x.center, x.scale) 
{
    x <- scale(x, center = F, scale = 1/x.scale)
    x <- scale(x, center = -x.center, scale = F)
    x
}
