unrollZGrid<- function( grid.list, ZGrid){
  if( is.null(ZGrid)){
    return(ZGrid)
  }
  if( is.list( ZGrid) ){
     if( any(grid.list[[1]] != ZGrid[[1]]) |any(grid.list[[2]] != ZGrid[[2]]) ){
         stop("grid list does not match grid for covariates")
       }  
# wipe out the x and y components of ZGrid because grid.list will be used
  ZGrid<- ZGrid$z
  }
# check dimensions
    Zdim<- dim( ZGrid)
      nx<- length( grid.list[[1]])
      ny<- length( grid.list[[2]])
      if( (Zdim[1] != nx) | (Zdim[2] != ny) ){
         stop( "Dimension of ZGrid does not match dimensions of location grid list.")
      }
# reshape as a matrix where rows index locations.
# Note that this works whether Zdim[3] exists or not! 
      return( matrix( c(ZGrid),  nrow= Zdim[1]*Zdim[2] ))
 }
