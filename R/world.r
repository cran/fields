"world" <-
function(ylim = c(-90, 90), xlim = NULL, add = FALSE, asp = 1, 
    xlab = "", ylab = "", xaxt = "n", yaxt = "n", eps = 0.1,col=1, 
    shift = FALSE, 
    fill=FALSE, col.water="white", col.land="darkgrey", ...) 
{
# load world dat set of land outlines 
# (should not reload if it is already loaded)
  
       data(world.dat)

# check some options   
   if( shift & fill){ 
         stop("filling in land not implemented with shift option")}


   if (shift) {
        ind1<- !is.na(world.dat$x)
        ind2<- (world.dat$x<0)
# shift coordinates
        world.dat$x[ind2&ind1] <- world.dat$x[ind2&ind1] + 360
# "pick up pen" at the new edges by inserting NA's
        world.dat$x[(world.dat$x <= eps | world.dat$x >= (360 - 
                 eps))&ind1] <- NA}

    if (is.null(xlim)) {
        if (shift) {
            xlim <- c(0, 360)
        }
        else {
            xlim <- c(-180, 180)
        }
    }

# create new plotting region if add is FALSE
    if (!add) {
        plot(world.dat, ylim = ylim, xlim = xlim, type = "n", 
            xaxt = xaxt, yaxt = yaxt, xlab = xlab, ylab = ylab, 
            asp = asp, ...)
    }

# decide whether to draw lines or fill in land masses and lakes.
    if( !fill){
       lines(world.dat, err = -1, col=col, ...) }
     else{
       world.color(col.water=col.water, col.land=col.land,...)}
 
 
    invisible()
}


world.color<- function( xlim=c(-180, 180), ylim=c(-90,90),col.water="white", 
                                col.land="darkgrey", ... )
{
# load world dat set of land outlines 
# (should not reload if it is already loaded)
              data(world.dat)

# logicals for land and lakes in case these need to 
# be modified.

land<- TRUE
lakes<- TRUE

# first add ocean color as background
# find the size of current region but don't let it get bigger than
#lon lat range

 rect(xlim[1],ylim[1], xlim[2],ylim[2], col=col.water, border=NA)

# find separate polygons of land masses and islands in world.dat 
# these are indicated by NA's

 ind<-  (1:length( world.dat$x))[ is.na( world.dat$x)]
 ind<- c( 1, ind)
 N<- length(ind) - 1 # number of distinct polygons in map

 lakes.id<- c(46, 53, 25,26, 28,27, 4,47,48,51,49, 50)
 land.id<- (1:N)[ -lakes.id]

# loop through polygons of land
 if(land){
  for(  k in land.id){
    tempi <- ind[k] : ind[k+1]
    polygon( 
      list( x= world.dat$x[tempi], y=world.dat$y[tempi] ), 
      col=col.land, border=col.land,...)}
 }
#
# add in interior of antarctica as a 
# rectangle below farthest south coastline

 ytemp<- world.dat$y[ ind[21]+1 ]
    rect( -180, -90, 180, ytemp,  col=col.land, border=col.land,...)

# loop through polygons of lakes do this second so 
# color over writes the land 
 if( lakes){
    for(  k in lakes.id){
       tempi <- ind[k] : ind[k+1]
       polygon(
       list( x= world.dat$x[tempi], y=world.dat$y[tempi] ),
       col=col.water, border=col.water,...)}
 }

}

