\name{plot.surface}
\alias{plot.surface}
\title{
  Plots a surface 
}
\description{
Plots a surface object in several different ways to give 3-d
information e.g. a contour plots, perspective plots.  
}
\usage{
plot.surface
(obj, main = NULL, type = "b", zlab = NULL, xlab = NULL,
    ylab = NULL, levels = NULL, zlim = NULL, graphics.reset = NULL,
   ...)
}
\arguments{
\item{obj}{
A surface object. At the minimum a list with components x,y and z 
in the same form as the input list for the standard contour and persp
functions. This can also be an object from predict.surface. 
}
\item{main}{
Title for plot. 
 
}
\item{type}{
type="p" for a perspective plot, type="c" for a contour plot  
and type="b" for both (default). type="I" for an image plot with a legend
strip (see image.plot). type="C" is the "I" option with contours lines
added. 
}
\item{zlab}{
z-axes label 
}
\item{xlab}{
x-axes label 
}
\item{ylab}{
y-axes labels 
}
\item{levels}{
Vector of levels to be passed to contour function. 
}
\item{graphics.reset}{
Reset to original graphics parameters after function plotting. 
Default is to reset if type ="b" but not for the single plot options.  
}
\item{\dots}{
Other graphical parameters that are passed along to either persp, contour
or image. 
}
}
\seealso{
surface, predict.surface, as.surface,   
}
\examples{
fit<- Tps( BD[,1:4], BD$lnya) # fit surface to data 

# make grid list 
grid<-list( KCl="x", MgCl2=mean(BD[,2]),KPO4="y", dNTP=mean(BD[,4]))

# evalute on a grid on two 
# variables holding two fixed at their mean levels 
out.p<-predict.surface(fit, grid, extrap=T)  

plot.surface(out.p, type="C") # surface plot  
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
