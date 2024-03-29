#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2022 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.edu,
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2
##END HEADER
##END HEADER


# test of vgram

suppressMessages(library(fields))
options(echo=FALSE)

data( ozone2)

y<- ozone2$y[16,]
x<- ozone2$lon.lat

vgram( x, y, lon.lat=TRUE)-> out

# compute "by hand"

outer( y, y ,"-")-> hold
hold<- .5*hold^2 
rdist.earth( x,x)-> hold2
col( hold2)> row( hold2)-> upper

hold<- hold[upper]
hold2<- hold2[upper]
 order( hold2)-> od
hold2<- hold2[od]
hold<- hold[od]
ind<- is.na(hold)
hold<- hold[!ind]
hold2<- hold2[!ind]

test.for.zero( hold, out$vgram, tag="vgram single time")


# multiple times including NAs at some times

y<- t(ozone2$y[16:18,])
x<- ozone2$lon.lat[,]

out<- vgram( x, y, lon.lat=TRUE)


N<- nrow( y)

hold<-  cbind(c(outer( y[,1], y[,1],"-")),
         c(outer( y[,2], y[,2],"-") ),
         c(outer(y[,3], y[,3],"-"))  )
hold<- .5*hold^2
hold<- rowMeans( hold, na.rm=TRUE)
hold<- matrix( hold, N,N)

rdist.earth( x,x)-> hold2

col( hold2)> row( hold2)-> upper
hold<- hold[upper]
hold2<- hold2[upper]

order( hold2)-> od
hold2<- hold2[od]
hold<- hold[od]

ind<- is.na(hold)
hold<- hold[!ind]
hold2<- hold2[!ind]

test.for.zero( hold, out$vgram, tag="vgram more than one time point")

# test covariogram versus correlogram
y<- ozone2$y[16,]
x<- ozone2$lon.lat

tau2 = var(y, na.rm=TRUE)
lookCov = vgram(x, y, lon.lat=TRUE, type="covariogram")
lookCor = vgram(x, y, lon.lat=TRUE, type="correlogram")

test.for.zero(lookCov$vgram*(1/tau2), lookCor$vgram, tag="correlogram versus covariogram")

# test cross-covariogram versus cross-correlogram

tau2 = var(y, na.rm=TRUE)
lookCov = crossCoVGram(x, x, y, y, lon.lat=TRUE, type="cross-covariogram")
lookCor = crossCoVGram(x, x, y, y, lon.lat=TRUE, type="cross-correlogram")

test.for.zero(lookCov$vgram*(1/tau2), lookCor$vgram, tag="correlogram versus covariogram")
