"US" <-
function (xlim = c(-124.7,  -67.1), ylim = c(25.2, 49.4), add = FALSE, ...) 
{
if( !exists("US.dat"))data(US.dat)
    if (!add) {
        plot(US.dat$x,US.dat$y, ylim = ylim, xlim = xlim, xlab = "", 
            ylab = "", type = "n", xaxt = "n", yaxt = "n", ...)
    }
    lines(US.dat$x,US.dat$y, err = -1, ...)
    invisible()
}
