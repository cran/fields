"world" <-
function (ylim = c(-90, 90), xlim = c(-180, 180), add = FALSE,
            asp = 1, xlab = "", ylab = "",  xaxt = "n", yaxt = "n", ...)
{
    if (!add) {
        plot(world.dat, ylim = ylim, xlim = xlim, type = "n",
             xaxt = xaxt, yaxt = yaxt, xlab = xlab, ylab = ylab,
             asp = asp, ...)
    }
    lines(world.dat, err = -1, ...)
    invisible()
}
