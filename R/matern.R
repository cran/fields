"matern" <-
function (t , scale = 1, range = 1,alpha=1/range, 
    smoothness = 0.5, nu= smoothness, phi=scale) 
{
#
# Matern covariance function transcribed from Stein's book page 31
# nu==smoothness, alpha ==  1/range
# GeoR parameters map to kappa==smoothness and phi == range 
# check for negative distances
if( any( t <0)) stop("distance argument must be nonnegative")

t<- t*alpha
#
# the hairy constant ...
# this is different from Stein to make this a correlation function when 
# scale =1
con<- (2^(nu - 1))*gamma(nu)
con<- 1/con
#
# call to  Bessel function from R base package
#
 phi*ifelse( t>0,
                 con*(t^nu) * besselK(t , nu ), 1)
}
