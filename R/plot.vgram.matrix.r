"plot.vgram.matrix" <-
function(obj,...){
# check if just radil distance has been used for vgram
collapse<- is.null( obj$ind)
#
if( !collapse){

nx<- max( obj$ind[,1])
ny<- max(  obj$ind[,2])
temp<-  matrix( NA,nrow=nx+1, ncol=ny+1)
temp[ obj$ind+1] <- obj$vgram
image.plot( 0:nx, 0:ny, temp, xlab="X", ylab="Y",...)
}
else( plot( obj$d, obj$vgram))

}
