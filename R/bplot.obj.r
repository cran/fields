"bplot.obj" <-
function (data, pos = NA, width=NULL, labels = NULL, srt =
NULL, 
    add = F, space = 0.25, sort.names = F, xlab = "", ylab = "", 
    label.cex = 1, xaxt = "n", outlier = T, horizontal=F,...) 
{
    cols <- length(data)
    range.data <- c(NA, NA)

    
    if (is.null(labels)) {
        labels <- names(data)
    }
#
# if pos is not passed then make them integers. 
#
    if (is.na(pos[1])) {
      pos <- 1:cols
#      
   if (sort.names) {
            pos <- order(labels)
        }

    }

    if (is.null(width)) {
        width <- min(diff(sort(pos))) * space
        if (cols == 1) 
            width <- space
    }
    if (length(width) == 1) 
        width <- rep(width, cols)
if (!add) {
# find overall range of boxplots
  for (k in 1:cols) {
        range.data <- range(c(range.data, data[[k]]$range), na.rm = T)
  }
         temp1<-range.data
         temp2 <-  range(c(pos - (0.5 * width)/space, pos + (0.5 *width)/space))

#
# setup plot for adding the box plots.
#

if( horizontal)
{plot(temp1, temp2,type = "n", yaxt = xaxt,xlab = "", ylab = "", ...)}
 else{plot(temp2, temp1,type = "n", xaxt = xaxt, xlab = "", ylab = "",...)}
}

#
# now add the box plots 
#
    
    for (i in 1:cols) {
# modify draw.bplot.obj to custmize how boxes are drawn 
       draw.bplot.obj(data[[i]], width[i], pos[i], outlier = outlier,
             horizontal=horizontal)
    }


# add labels if they exist
 if (label.cex > 0) {
        if (is.null(srt)) {
# keep the labels from running into each other. 
#rotate labels if there are more than 7

            if (length(labels) > 7) {
                srt <- 90
            }
            else {
                srt <- 0
            }
        }
      if( horizontal){axis.loc<- 2}
      else{axis.loc<- 1}

         axis(axis.loc, pos, labels, tick = F, srt = srt, adj = .5, 
            cex = label.cex)

}
    invisible()
}
