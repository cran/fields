\name{pushpin}
\alias{pushpin}
%- Also NEED an '\alias' for EACH other topic documented here.
 \title{ Adds a "push pin" to an existing  3-d plot} 
 \description{ Adds  to an existing 3-d perspective plot 
a push pin to locate a  specific point.}
}
\usage{
pushpin( x,y,z,p.out, height=.05,col="black",text=NULL,adj=-.1,cex=1.0,...)

}

\arguments{
  \item{x}{x location}
  \item{y}{y location}
  \item{z}{z location}
  \item{p.out} {Projection information returned by persp}   
  \item{height}{Height of pin in device coordinates (default is about 5\%). }
  \item{col}{Color of pin head}
  \item{text}{Optional text to go next to pin head}
  \item{adj}{Position of text relative to pin head.}
  \item{cex}{Character size for in head and/or text}
  \item{\dots}{Additional graphics arguments that are passed to the text 
                     function.}
}

\details{
See the help(text) for the conventions on
the \code{adj} argument and other options for placing text.  

 }
\author{Doug Nychka}

\seealso{drape.plot,persp}

\examples{
# Dr. R's favorite New  Zealand Volcano!
     data( volcano)
     M<- nrow( volcano)
     N<- ncol( volcano)
     x<- seq( 0,1,,M)
     y<- seq( 0,1,,N)

     drape.plot( x,y,volcano, col=terrain.colors(128))-> pm 

     max( volcano)-> zsummit
     xsummit<- x[ row( volcano)[volcano==zsummit]]
     ysummit<- y[ col( volcano)[volcano==zsummit]]

pushpin( xsummit,ysummit,zsummit,pm, text="Summit")


}

\keyword{ hplot }% at least one, from doc/KEYWORDS
