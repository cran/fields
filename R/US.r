# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"US" <-
function (xlim = c(-124.7, -67.1), ylim = c(25.2, 49.4), add = FALSE, 
    shift = FALSE, ...) 
{
    if (!exists("US.dat")) 
        data(US.dat)
    if (shift) {

    ind1<- !is.na(US.dat$x)
    ind2<- US.dat$x < 0
    US.dat$x[ind2&ind1] <- US.dat$x[ind2&ind1] + 360
    xlim<- c(-124.7, -67.1) + 360
    }
    if (!add) {
        plot(US.dat$x, US.dat$y, ylim = ylim, xlim = xlim, xlab = "", 
            ylab = "", type = "n", xaxt = "n", yaxt = "n", ...)
    }
    lines(US.dat$x, US.dat$y, err = -1, ...)
    invisible()
}

