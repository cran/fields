%#
%# fields  is a package for analysis of spatial data written for
%# the R software environment.
%# Copyright (C) 2026 Colorado School of Mines
%# 1500 Illinois St., Golden, CO 80401
%# Contact: Douglas Nychka,  douglasnychka@gmail.com,
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

\name{BD}
\alias{BD}
\title{
 Effect of buffer compositions on DNA strand displacement
amplification. 
}
\description{
A 4-d regression data set
with replication. This is a useful test data set for function
fitting methods.
\code{BD} as a data frame has 89 rows and 5 columns. There are 89 runs with
four buffer components (KCl, MgCl2, KPO4, dNTP) systematically varied in
a space-filling design. The response is the DNA amplification rate.
}
\usage{ data(BD) }
\format{
This data frame contains the following columns:

\describe{
\item{KCl}{
Buffer component.
}
\item{MgCl2}{
Buffer component.
}
\item{KPO4}{
Buffer component.
}
\item{dNTP}{
Buffer component, deoxyribonucleotides.
}
\item{lnya}{
Exponential amplification rate on a log scale, i.e. the actual amplification
rate.
}
}
}

\source{
Thanks to Perry Haaland and Michael O'Connell for these data.
at Becton Dickinson Research Center,
Research Triangle Park, NC
}
\seealso{
Tps, spatialProcess 
}
\examples{
# fitting a DNA strand
# displacement amplification  surface to various buffer compositions
data(BD)
fit<- Tps(BD[,1:4],BD$lnya,scale.type="range")
surface(fit)  # plots fitted surface and contours
}
\keyword{datasets}
% docclass is function
% Converted by Sd2Rd version 1.21.
