
color.scale<- function( z, col=tim.colors(256), zlim =NULL, transparent.color="white",
                       eps= 1e-8){
#  
# converts real values to a color scale of NC values. 
# role of eps is to prevent values exactly at the end of the range from being
# missed
  if( is.null( zlim)){ zlim <- range(z, na.rm=TRUE)}
    z[(z < zlim[1]) | (z > zlim[2])] <- NA
   NC <- length(col)
   breaks<- seq(zlim[1]*(1-eps), zlim[2] * (1 + eps) ,, NC+1)
# the magic of R ... 
    icolor<-  cut(c(z),breaks)@.Data
# returned values is a vector of character hex strings encoding the colors.
  ifelse(is.na(icolor), transparent.color, col[icolor])
}
  
