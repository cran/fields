# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

Krig.null.function<- function( x, Z=NULL,drop.Z=FALSE, m){

# default function to create matrix for fixed part of model
#  x, Z, and drop.Z are required
#  Note that the degree of the polynomial is by convention (m-1)
# returned matrix must have the columns from Z last!
# 
  if( is.null( Z)| drop.Z){
    
     return(fields.mkpoly( x, m=m))}
  else{
     return(cbind(fields.mkpoly( x, m=m),Z)) }

}
