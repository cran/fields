"world" <-
function (ylim = c(-90, 90), xlim = c(-180, 180), add = F, ...) 
{
    if (!add) {
        plot(world.dat, ylim = ylim, xlim = xlim, xlab = "", 
            ylab = "", type = "n",  xaxt = "n", yaxt = "n", 
            ...)
    }
    lines(world.dat, err = -1, ...)
    invisible()
}
