% fields, Tools for spatial data
% Copyright 2004-2011, Institute for Mathematics Applied Geosciences
% University Corporation for Atmospheric Research
% Licensed under the GPL -- www.gpl.org/licenses/gpl.html

\name{bplot}
\alias{bplot}
\title{
  boxplot  
}
\description{
Plots boxplots of several groups of data
and allows for placement at different horizontal or vertical positions or 
colors. It
is also flexible in the input object, accepting either a list or matrix.
}
\usage{
bplot(x, by, pos=NULL, at = pos, add = FALSE, boxwex =
                 0.8,xlim=NULL, ...)
}
\arguments{
 \item{x}{
Vector, matrix,  list or data frame. A vector may be divided according
to the by argument. Matrices and data frames are separated by columns and
lists by components.   
}
\item{by}{
If x is a vector, an optional vector (either character or numerical)
specifying the categories to divide x into separate data sets. Boxplots are then made 
for each group.
}
 \item{pos}{
  The boxplots will be plotted vertically (horizontally) and pos gives the x (y) locations for their centers. If omitted the boxes are equally spaced at integer  
  values.  This is the same as \code{at} in the \code{boxplot} function 
  }
\item{at}{Same as \code{pos} this is the name for this argument in the standard \code{boxplot} function.}

\item{add}{
   If true, do not create a new plots just add the boxplots to a current
   plot. Note that the pos argument may be useful in this case and should 
   be in the user coordinates of the parent plot.}
\item{boxwex}{A boxplot argument to control the width of the boxplot.
It behaves a little different than as an argumetn passed directly to \code{boxplot}.
To make this a general function it is useful to scale this according to size of positions. Within bplot this happens as \code{boxwex<- boxwex* min(diff( sort( at)))}.
and then the scaled version of \code{boxwex} is now passed to \code{boxplot}.}
\item{xlim}{ Same as the usual argument used in plotting. The plotting limits for the 
x axis. }
\item{\dots}{
Other arguments to be passed to the boxplot function some handy favorites are:
\code{names}
  Labels for each boxplot.
\code{horizontal}If TRUE draw boxplots horizontally the default is false, produce
    vertical box plots.
\code{lwd}Width(s) of lines in box plots. 
\code{col}Color(s) of bplots. See \code{colors()} for some choices.}   
}

\details{
 This function was created as a complement to the usual S function for
boxplots. The current function makes it possible to put the boxplots
at unequal x or y positions in a rational using the \code{at} or
\code{pos} arguments. This is useful for visually grouping a large set
of boxplots into several groups. Also placement of the boxplots with
respect to the axis can add information to the plot.  Another aspect
is the emphasis on data structures for groups of data. One useful
feature is the by option to break up the x vector into distinct
groups.


The older \code{bplot} function drew the boxplots from scratch and if
one needs to do this refer to the old functions: \code{
describe.bplot, draw.bplot.obj, bplot.xy, bplot.obj}

Finally to bin data into groups based on a continuous variable and to
make bplots of each group see \code{bplot.xy}.
}
\seealso{ bplot.xy }
\examples{
#
set.seed(123)
temp<- matrix( rnorm(12*8), ncol=12)
pos<- c(1:6,9:14)
bplot(temp)
#
par(las=1)
bplot( temp, pos=pos, names=paste( "Data",1:12, sep=""), horizontal=TRUE)

#
# Xmas boxplots in red and green 
bplot( temp, pos=pos,  col=c("red", "green"))
# add an axis on top
axis( 3)
}
\keyword{hplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
