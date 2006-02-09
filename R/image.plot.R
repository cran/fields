"image.plot" <-
function (..., add = FALSE, nlevel = 64, legend.shrink = 0.9, 
    legend.width = 1.2, legend.mar = NULL, graphics.reset = FALSE, 
    horizontal = FALSE, bigplot = NULL, smallplot = NULL, legend.only = FALSE, 
    col = tim.colors(nlevel)) 
{
    old.par <- par(no.readonly = TRUE)

#  figure out zlim from passed arguments

    info <- image.plot.info(...)

    if (add) {
        big.plot <- old.par$plt
    }
  
    if (legend.only) {
        graphics.reset <- TRUE
    }
  
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1,5.1)
    }
#
# figure out how to divide up the plotting real estate. 
# 
    temp <- image.plot.plt(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, legend.mar = legend.mar, 
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
#
# bigplot are plotting region coordinates for image
# smallplot are plotting coordinates for legend 
 
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot

#
# draw the image in bigplot, just call the R base function 
#
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        image(..., add = add, col = col)
        big.par <- par(no.readonly = TRUE)
    }

##
## check dimensions of smallplot

    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }

# Following code draws the legend using the image function 
# and a one column image. 
    
#
# OLD code misaligns     iy <- seq(info$zlim[1], info$zlim[2], , nlevel)
# This is corrected -- thanks to S. Woodhead
#
    ix <- 1
    minz <- info$zlim[1]
    maxz <- info$zlim[2]
    binwidth <- (maxz - minz) / nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints

    iz <- matrix(iy, nrow = 1, ncol = length(iy))

    breaks<- list(...)$breaks

# draw either horizontal or vertical legends. 
# either using breaks or not 

    if (!horizontal) {
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)

        if( is.null( breaks)){
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col )}
        else{
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                 ylab = "", col = col, breaks=breaks )}

        axis(4, mgp = c(3, 1, 0), las = 2)
        box()
    }
    else {
        par(new = TRUE, pty = "m", plt = smallplot, err = -1)
      
        if( is.null( breaks)){
            image(iy, ix, t(iz), yaxt = "n", xlab = "", ylab = "", 
                col = col)}
        else{
            image(iy, ix, t(iz), yaxt = "n", xlab = "", ylab = "",
                col = col,breaks=breaks)}

        box()
    }

# clean up graphics device settings


    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}

