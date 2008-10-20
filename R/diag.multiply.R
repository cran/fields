# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html


# create a special matrix multiply for diagonal
# matrices. Diagonal matrix assumed to be just a vector.
# see tests directory for tests
# NOTE: this is not a symmetric operation:
#  when a left vector is given it is a diagonal matrix
#  when a right vector is given it is a vector. 
#

setGeneric("%d*%",function(x,y,...)standardGeneric("%d*%"))

setMethod("%d*%",signature(x="matrix",y="matrix"),
function(x,y){x%*%y} )

setMethod("%d*%",signature(x="matrix",y="numeric"),
function(x,y){
  if( length(y)!= ncol(x)){
    stop("non-conformable arguments in %d*% ")}
    matrix(t(y*t(x)), ncol=ncol(x), nrow=nrow(x))} )

setMethod("%d*%",signature(x="numeric",y="matrix"),
function(x,y){
  if( length(x)!= nrow(y)){
    stop("non-conformable arguments in %d*%")}
    matrix(x*y, ncol=ncol(y), nrow=nrow(y))} )

setMethod("%d*%",signature(x="numeric",y="numeric"),
function(x,y){cbind(x*y)} )
