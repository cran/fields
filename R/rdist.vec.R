# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
rdist.vec = function(x1, x2) {
  #make sure inputs are matrices
  if (!is.matrix(x1)) {
    x1 <- as.matrix(x1)
  }
  if(!is.matrix(x2)) {
    x2 <- as.matrix(x2)
  }
  
  #return distances
  sqrt(rowSums((x1 - x2)^2))
}