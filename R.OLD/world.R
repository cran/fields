# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"world" <- function( ...){
     map("world", ...)
    invisible()
}
world.color <- function(...){
 cat("world.color has been depreciated. Please use fill options in
the world/map function.", fill=TRUE)
}
world.land <- function( ...){
 cat("world.land has been depreciated. Please use fill options in
the world/map function.", fill=TRUE)
}

in.land.grid<- function(...)
{
 cat("world.land has been depreciated. Please refer to fields 6.7.1 or earlier to acces this function.", fill=TRUE)
}


