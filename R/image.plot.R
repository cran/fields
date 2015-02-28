# fields, Tools for spatial data
# Copyright 2004-2013, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"image.plot" <- function(..., add = FALSE,
    breaks= NULL, nlevel = 64, col = NULL,  
    horizontal = FALSE, legend.shrink = 0.9, legend.width = 1.2, 
    legend.mar = ifelse(horizontal, 3.1, 5.1), legend.lab = NULL,
    legend.line= 2,                    
    graphics.reset = FALSE, bigplot = NULL, smallplot = NULL, 
    legend.only = FALSE,  lab.breaks = NULL, 
    axis.args = NULL, legend.args = NULL, legend.cex=1.0, midpoint = FALSE, border = NA, 
    lwd = 1, verbose=FALSE) {
    # Thanks to S. Koehler and  S. Woodhead
    # for comments on making this a better function
    #
    # save current graphics settings
    old.par <- par(no.readonly = TRUE)
    # set defaults for color scale 
    # note this works differently than the image function.
    if( is.null(col))  {
    	col<-  tim.colors(nlevel)}
    	else{
    		nlevel<- length( col)
    		}
    #  figure out zlim from passed arguments
    #  also set the breaks for colors if they have not been passed, 
    info <- imagePlotInfo(..., breaks=breaks, nlevel=nlevel)
    # breaks have been computed if not passed in the call
    breaks<- info$breaks
    if( verbose){
    	print(info)
    }
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    # figure out how to divide up the plotting real estate 
    temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
        legend.width = legend.width, legend.mar = legend.mar, 
        horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    # bigplot has plotting region coordinates for image
    # smallplot has plotting coordinates for legend strip
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    # draw the image in bigplot, just call the R base function
    # or poly.image for polygonal cells
    # note the logical switch
    # for poly.grid is parsed out of call from image.plot.info
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        if (!info$poly.grid) {
            image(..., breaks=breaks, add = add, col = col)
        }
        else {
            poly.image(..., add = add, col = col, midpoint = midpoint, 
                border = border, lwd.poly = lwd)
        }
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
    # What might be confusing is the values of the "image" are the same 
    # as the locations on the legend axis.
    # Moreover the image values are in the middle of each breakpoint category
    # thanks to Tobias Nanu Frechen and Matthew Flickinger 
    # for sorting out some problems with the breaks position in the legend.
    ix <- 1:2
    iy<- breaks
    nBreaks<- length( breaks)
    midpoints<- (breaks[1:(nBreaks-1)] +  breaks[2:nBreaks] )/2
    iz <- matrix(midpoints, nrow = 1, ncol = length(midpoints)) 
    if( verbose){print(breaks)
    	print( midpoints)
    	print( ix)
    	print( iy)
    	print( iz)
    	print( col)}  
        #
    # next par call sets up a new plotting region just for the legend strip
    # at the smallplot coordinates
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    # draw color scales the two  cases are horizontal/vertical 
    # add a label if this is passed.
    if (!horizontal) {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks=breaks)
    }
    else {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                ylab = "", col = col, breaks=breaks)
    }
    # create the argument list to draw the axis
    #  this avoids 4 separate calls to axis and allows passing extra
    # arguments.
    if (!is.null(lab.breaks)) {
        # axis with labels at break points
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
            at = breaks, labels = lab.breaks), axis.args)
    }
    else {
        # If lab.breaks is not specified ( with or without breaks), pretty
        # tick mark locations and labels are computed internally,
        # or as specified in axis.args at the function call
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
            axis.args)
    }
    #
    # now add the axis to the legend strip.
    # notice how all the information is in the list axis.args
    do.call("axis", axis.args)
    # add a box around legend strip
    box()
    #
    # add a label to the axis if information has been  supplied
    # using the mtext function. The arguments to mtext are
    # passed as a list like the drill for axis (see above)
    #
    if (!is.null(legend.lab)) {
        legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
            1, 4), line = legend.line, cex=legend.cex)
        #                    just guessing at a good default for line argument!
    }
    # add the label using mtext function
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
    #
    # clean up graphics device settings
    # reset to larger plot region with right user coordinates.
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
        # Suggestion from Karline Soetaert <Karline.Soetaert@nioz.nl>
        # this is to reset margins to be based on the mar arguments
  #      par(mar = par("mar"))  or
  #      par(mar = big.par$mar)
  # unfortunately this causes problems by allowing plotting outside of the
  # original plot region.
        invisible()
    }
}
