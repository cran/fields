"world" <-
function (ylim = c(-90,90), xlim = NULL, add = FALSE, 
    asp = 1, xlab = "", ylab = "", xaxt = "n", yaxt = "n", 
eps=.1,shift=FALSE,...) 
{

if( shift){
# this changes
#the range of lon from (-180, 180) to (0,360)
world.dat$x[world.dat$x<0] <- world.dat$x[world.dat$x<0] +360
#
# find where there are lines that cross lon=0
# and add a break using NA

world.dat$x[world.dat$x<=eps|world.dat$x>=(360-eps)]<-NA

} 

if( is.null(xlim) ) {
   if( shift)
            {xlim<- c(0,360)}
         else
            {xlim<- c(-180,180)}
    }

    if (!add) {
        plot(world.dat, ylim = ylim, xlim = xlim, type = "n", 
            xaxt = xaxt, yaxt = yaxt, xlab = xlab, ylab = ylab, 
            asp = asp, ...)
    }
    lines(world.dat, err = -1, ...)
    invisible()
}
