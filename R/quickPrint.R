quickPrint<- function( obj, max.values=10){
# this simple function will only print an object if
# it is not too big	
	if( is.list( obj)){
	 sizeObj<- length( unlist( obj))			
	}
	else{
	 sizeObj<-length( c( obj))
	}
	if( sizeObj<= max.values){
		print(obj)
	}
	else{
		cat("Object size is more than", max.values,"items (", sizeObj, "total)",  fill=TRUE)
	}
}
