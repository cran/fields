"image.plot" <-
function (..., add = F, nlevel = 32, legend.shrink = 0.9, legend.width = 0.04, 
    graphics.reset = F, horizontal = F, offset = 2 * legend.width, 
    bigplot = NULL, smallplot = NULL, legend.only = F, col = topo.colors(64)) 
{
    old.par <- par(no.readonly = T)
    info <- image.plot.info(...)
    if (add) 
        big.plot <- old.par$plt
    if (legend.only) 
        graphics.reset <- T
    temp <- image.plot.plt(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, horizontal = horizontal, 
        offset = offset, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        image(..., add = add, col = col)
        big.par <- par(no.readonly = T)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    temp <- list(...)
    iy <- seq(info$zlim[1], info$zlim[2], , nlevel)
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    ix <- 1
    if (!horizontal) {
        par(new = T, pty = "m", plt = smallplot, err = -1)
        image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
            ylab = "", col = col)
        axis(4, mgp = c(3, 0, 0))
    }
    else {
        par(new = T, pty = "m", plt = smallplot, err = -1)
        image(iy, ix, t(iz), yaxt = "n", xlab = "", ylab = "", 
            col = col)
    }
mfg.save<- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
par( mfg=mfg.save, new=F)
        invisible()
    }
    else {
        par(big.par)
par( mfg=mfg.save, new=F)
        invisible()
    }
}
