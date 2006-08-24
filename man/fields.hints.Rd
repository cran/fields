\name{fields.hints}
\alias{fields.hints}
\alias{fields.style}

\title{
fields - graphics hints  
}

\description{
 Here are some technical hints for assembling multiple plots with common
legends or axes and setting the graphics parameters to make more
readable figures . These use the standard graphics environment. 
}
\usage{
fields.style()
}
\details{
\code{fields.style} is a simple 
function to enlarge the characters in a plot and set the colors. List this out to modify the choices.

}

\examples{


##Two examples of concentrating a panel of plots together
## to conserve the white space. 
## see also the example in image.plot using split.screen. 
## The basic trick is to use the oma option to reserve some space around the 
## plots.  Then unset the outer margins to use that room. 

library( fields)

# some hokey image data
x<- 1:20
y<- 1:15
z<-  outer( x,y,"+")
zr<- range( c(z))

# add common legend to 3X2 panel 

par( oma=c(4,0,0,0))
set.panel( 3,2)
par( mar=c(1,1,0,0))

# squish plots together with just 1 space between
for ( k in 1:6){
image( x,y,z, axes=FALSE, xlab="", ylab="", zlim=zr)
}

par( oma=c(0,0,0,0))
image.plot( zlim=zr, legend.only=TRUE, horizontal=TRUE, legend.mar=5)

# you may have to play around with legend.mar and the oma settings to
# get enough space.



##
### also add some axes on the sides. in a lattice style 
## note oma adds some more room at bottom. 

par( oma=c(8,6,1,1))
set.panel( 3,2)
par( mar=c(1,1,0,0))
##
for ( k in 1:6){
 image( x,y,z, axes=FALSE, xlab="", ylab="", zlim=zr)
 box() # box around figure

# maybe draw an x axis 
  if( k  \%in\% c(5,6) ){
  axis( 1, cex.axis=1.5)
  mtext( line=4, side=1, "Xstuff")}

# maybe draw a y axis
  if( k  \%in\% c(1,3,5) ){
  axis( 2, cex.axis=1.5)
  mtext( line=4, side=2, "Ystuff")}
}

# same trick of adding a legend strip. 
par( oma=c(0,0,0,0))
image.plot( zlim=zr, legend.only=TRUE, horizontal=TRUE, legend.mar=5)



####
# show colors
## the factory colors:

clab<- colors()
n<- length( clab)
N<- ceiling( sqrt(n) )
M<- N
temp<- rep( NA,M*N)
temp[1:n] <- 1:n
z<- matrix(temp, M,N)

image(seq(.5,M+.5,,M+1), seq(.5,N+.5,,N+1)
       , z,  col=clab, axes=FALSE, xlab="", ylab="")


# function to locate colors
#
show.color<- function(ind=NA){
 c( mar=c(0,0,3,0) )

 clab<- colors()
 n<- length( clab)
 N<- ceiling( sqrt(n) )
 M<- N
 temp<- rep( NA,M*N)
 temp[1:n] <- 1:n
 z<- matrix(temp, M,N)

 image(seq(.5,M+.5,,M+1), seq(.5,N+.5,,N+1)
       , z,  col=clab, axes=FALSE, xlab="", ylab="")

 if( is.na(ind)){
  cat("Use mouse to identify color",fill=TRUE)
  locator( 1)-> loc; i<- round( loc$x);j<- round( loc$y);
  ind<- z[i,j] }

 cat("ID ", ind," name ",clab[ind],fill=TRUE)
 cat( "RGB", col2rgb(clab[ind])/256, fill=TRUE)

 points( i,j,  col= clab[ind], cex=4, pch="O")
 points( i,j, pch="+", col= "black", cex=1)

 mtext( side=3, text=clab[ind], col=clab[ind], line=1, cex=2)}



}


\keyword{hplot}
