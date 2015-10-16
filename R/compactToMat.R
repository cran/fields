# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
compactToMat = function(compactMat, diagVal=0, lower.tri=FALSE, upper.tri=TRUE) {
  #compactMat: a symmetric matrix stored as a vector containing elements for the upper triangle
  #portion of the true matrix
  #diagVal: a number to put in the diagonal entries of the output matrix
  #lower.tri: if TRUE, fills in lower tringular portion of the matrix
  #upper.tri: if TRUE, fills in upper tringular portion of the matrix
  
  if(class(compactMat) == 'dist') {
    n = attr(compactMat, "Size")
  } else { # (n^2 - n)/2 = length(compactMat)
    stop("input matrix is not compact or is not of class \"dist\"")
    
    #or if class is not dist but input matrix is still compact, use:
    #n = (1 + sqrt(1 + 8*length(compactMat)))/2
  }
  
  return(.Call("compactToMatC", as.double(compactMat), length(compactMat), as.integer(n), as.double(diagVal), 
               as.integer(lower.tri), as.integer(upper.tri)))
}
