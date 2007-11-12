"crop.image" <-function(obj, loc=NULL,...){

 if( is.null(loc)){
    image.plot( obj,...)
    get.rectangle()-> loc}

# coerce to midpoints
  m<- nrow( obj$z)
  n<- ncol( obj$z)
  nx<- length( obj$x)
  ny<- length( obj$y)

  if( nx != m) { obj$x<- (obj$x[1:m] + obj$x[2:(m+1)])/2 }
  if( ny != n) { obj$y<- (obj$y[1:n] + obj$x[2:(n+1)])/2 }


# coerce loc to x,y list format if matrix  or data frame
 if( is.matrix( loc)|is.data.frame(loc) ){
       if( ncol(loc)!=2){  stop("loc must have two columns 
          (for x and y coordinates )")}
       loc<- list( x= loc[,1], y=loc[,2] ) }

 x<-obj$x
 y<- obj$y

 

 N<- length(x)
 xr<- range( loc$x)
 xtest<- range( x)
 if( xr[1] < xtest[1] | xr[2] > xtest[2]){ 
        stop("cropping outside ranges of x values")}

 x1<-  max( (1:N)[ xr[1]>=x])
 x2<-  min( (1:N)[xr[2]<=x])

 N<- length( y)
 yr<- range( loc$y)
 ytest<- range( y)
 if( yr[1] < ytest[1] | yr[2] > ytest[2]){
       stop("cropping outside ranges of y values")}

 y1<-  max( (1:N)[ yr[1]>=y])
 y2<-  min((1:N)[yr[2]<=y])

 list( x=obj$x[x1:x2], y=obj$y[y1:y2], z=obj$z[x1:x2,y1:y2])

}

"get.rectangle" <-
function(){
locator( 2, type="p", pch="+")-> temp
rect( temp$x[1], temp$y[1], temp$x[2], temp$y[2])
temp
}


"half.image" <-
function(obj){

# coerce to list if a matrix 
if( is.matrix(obj)){
obj<- list( x= 1:nrow( obj), y= 1:ncol(obj), z= obj)
}
 
M<- length( obj$x)
N<- length(obj$y)
M2<- trunc( M/2)
N2<- trunc( N/2)
z<- matrix( NA, nrow= M2, ncol=N2)
ix<- (1:M2)*2
iy<- (1:N2)*2 
x2<- (obj$x[ ix-1] + obj$x[ ix])/2 
y2<- (obj$y[ iy-1] + obj$y[ iy])/2 
return( 
  list( x=x2, y=y2, z=
     (obj$z[ix-1,iy] + obj$z[ix-1, iy-1] + obj$z[ix,iy-1] + obj$z[ix, iy])/4  )
)

}

pushpin<- function( x,y,z,p.out, height=.05,
                col="black",text=NULL,adj=-.1, cex=1.0,...){

trans3d( x,y,z, p.out)-> Sxy1

#trans3d( x,y,z, p.out)-> Sxy2

Sxy2<- Sxy1
par()$usr-> hold
Sxy2$y <- (hold[4]-hold[3])*height + Sxy2$y

segments( Sxy1$x, Sxy1$y, Sxy2$x, Sxy2$y, col="black")
points( Sxy2, col=col, pch=19, cex=cex)
if( !is.null(text))
{ text( Sxy2$x, Sxy2$y, label=text, adj=adj, cex=cex,...)}
}

"two.colors" <-
function (n=256, 
 start="darkgreen" , 
end="red", middle="white")
{
n1<- n/2
n2<- n- n1
col2rgb( end)-> col2 
col2rgb( start)-> col1 
col2rgb(middle)-> mid.col
e1<- seq( 1,0,,n1)
e2<- seq( 0,1,,n2)

temp<-rbind( e1* matrix( col1, nrow=n1, ncol=3,byrow=TRUE) 
       + (1-e1)* matrix( mid.col, nrow=n1, ncol=3,byrow=TRUE),

             e2* matrix( col2, nrow=n1, ncol=3,byrow=TRUE) 
       + (1-e2)* matrix( mid.col, nrow=n1, ncol=3,byrow=TRUE))

temp<- temp/256
rgb( temp[,1], temp[,2], temp[,3])

}

