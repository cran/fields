"image.plot" <-
function (..., add = FALSE, nlevel = 64, legend.shrink = 0.9, 
    legend.width = 0.05, graphics.reset = FALSE, horizontal = FALSE, 
    offset = 2 * legend.width, bigplot = NULL, smallplot = NULL, 
    legend.only = FALSE, col = topo.colors(nlevel)) 
{
    old.par <- par(no.readonly = TRUE)
    info <- image.plot.info(...)
    if (add) 
        big.plot <- old.par$plt
    if (legend.only) 
        graphics.reset <- TRUE
#
# compute what the layout should be. 
#
    temp <- image.plot.plt(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, horizontal = horizontal, 
        offset = offset, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {

        if (!add) {
            par(plt = bigplot)
        }
#
# the call to the good 'ole image function
#
        image(..., add = add, col = col)
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
#
# having no drawn the image in the right size plotting region 
# switch to what is left over and add the legend
# note the legend strip is actually just another skiiny image plot!
#
    temp <- list(...)
    iy <- seq(info$zlim[1], info$zlim[2], , nlevel)
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    ix <- 1

    if (!horizontal) {
# this is a vertical legend
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
        axis(4, mgp = c(3, 0, 0))
    }
    else {
# this is a horizontal legend
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
        image(iy, ix, t(iz), yaxt = "n", xlab = "", ylab = "", 
            col = col)
    }
#
# reset some of the graphics parameters. 
#
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
# this resets plotting region to the larger image plot
#      I am not sure why the the second line is  needed but this 
#      seems to set the right plotting limits  and make sure clipping 
#      occurs...
        par(big.par)
        par( plt=big.par$plt, xpd=FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}
