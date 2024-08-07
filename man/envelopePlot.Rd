%#
%# fields  is a package for analysis of spatial data written for
%# the R software environment.
%# Copyright (C) 2024 Colorado School of Mines
%# 1500 Illinois St., Golden, CO 80401
%# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
%#
%# This program is free software; you can redistribute it and/or modify
%# it under the terms of the GNU General Public License as published by
%# the Free Software Foundation; either version 2 of the License, or
%# (at your option) any later version.
%# This program is distributed in the hope that it will be useful,
%# but WITHOUT ANY WARRANTY; without even the implied warranty of
%# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%# GNU General Public License for more details.
%#
%# You should have received a copy of the GNU General Public License
%# along with the R software environment if not, write to the Free Software
%# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
%# or see http://www.r-project.org/Licenses/GPL-2
%##END HEADER
%##END HEADER
\name{envelopePlot}
\alias{envelopePlot}
\title{ Add a shaded the region between two functions to an existing plot}

\description{
This function shades the region vertically between two functions, specified as pairs of x and y vectors,
and draws the functions in a darker shade. More formally, it shades all points (x,y) such that f1(x) < y
< f2(x) or f2(x) < y < f1(x). When both functions have the same group of x values, 
the x values only need to be set once but y2 needs to be passed in by name. 
If the two functions intersect, the vertical space between the functions will be shaded on both sides, as
implied in the definition above.
}

\usage{
envelopePlot(x1, y1, x2 = x1, y2,
                         col ="thistle1" , lineCol = "thistle3", ...)
}

\arguments{
\item{x1}{The x coordinates for the first function (or possibly both functions).}
\item{y1}{The y coordinates for the first function.}
\item{x2}{The x coordinates for the second function.}
\item{y2}{The y coordinates for the second function.}
\item{col}{The color to make the filling between the functions.}
\item{lineCol}{The color to make the lines representing the functions.}
\item{\dots}{Additional arguments to the base R function \code{polygon} }
}

\author{
Matt Iverson
}

\examples{
x <- seq(0, 2*pi,, 100)
y1 <- cos(x)
y2 <- sin(x)
plot(x, y1, type="l")
envelopePlot(x, y1, y2=y2)

x1 <- c(0, 0.5, 1)
y1 <- c(0, 2, 1)
x2 <- c(0, 1)
y2 <- c(-1, 0)
plot(x1, y1, type="l", ylim=c(-1, 2))
envelopePlot(x1, y1, x2, y2)
}
