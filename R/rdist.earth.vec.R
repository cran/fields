# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
rdist.earth.vec = function(x1, x2, miles=TRUE, R=NULL) {
  
  #set default radius
  if(is.null(R)) {
    if(miles)
      R = 3963.34
    else
      R = 6378.388
  }
  
  #convert lon/lat to radians
  x1 = x1 * (pi/180)
  x2 = x2 * (pi/180)
  
  #calculate distances using Haversine method
  lonDist2 = (x2[,1] - x1[,1]) * (1/2)
  latDist2 = (x2[,2] - x1[,2]) * (1/2)
  a = sin(latDist2) * sin(latDist2) + cos(x1[, 2]) * cos(x2[, 2]) * sin(lonDist2) * sin(lonDist2)
  return(2 * atan2(sqrt(a), sqrt(1 - a)) * R)
}