# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"summary.ncdf" <-
function( object,...){

for ( i in 1:object$nvars ){

   vname = object$var[[i]]$name    # variable name
   ndims = object$var[[i]]$ndims   # number of dimensions for this variable

   dimstring = paste(vname,'( variable ',as.character(i),') has shape')
   for (j in 1:ndims) {
      dimstring <- paste(dimstring, 
        as.character(object$var[[i]]$dim[[j]]$len))
   }

   cat(dimstring, fill=TRUE)
}
}
