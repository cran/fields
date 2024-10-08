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

\name{xline}
\alias{xline}
\title{
  Draw a vertical line 
}
\description{
Adds  vertical lines in the plot region. 
}
\usage{
xline(x, ...)
}
\arguments{
\item{x}{
Values on x axis specifying location of vertical lines. 
}
\item{\dots}{
Any ploting options for abline. 
}
}
\seealso{
 \link{yline}, \link{abline} 
}
\examples{

plot( 1:10)
xline( 6.5, col=2)
 
world( col=3) 
yline( seq( -80,80,10),col=4, lty=2)
xline( seq( -180,180,10),col=4,lty=2)
yline( 0, lwd=2, col=4)
}
\keyword{aplot}
% docclass is function
% Converted by Sd2Rd version 1.21.
