"matern.earth.cov" <-
function (x1, x2, theta = 1, smoothness = 0.5, 
             scale = 1, miles = TRUE, R = NULL) {

    if (!is.matrix(x1)) 
        stop( "x1 must a two column matrix of lat/lons")
    if (missing(x2)) 
        x2 <- x1
    if (!is.matrix(x2)) 
        stop( "x1 must a two column matrix of lat/lons")
    if (length(theta) != 1) 
        stop(" theta must be a scalar")
      
      if (is.null(R)) {
         if (miles) 
            R <- 3963.34
        else R <- 6378.388
      }

# separation in radians
    temp<-  rdist.earth( x1,x2, R=R, miles=miles)/R

# strange transformation to guarentee positive definiteness
     temp<- 2*R* sin( temp/2 )

# now apply the "usual" Matern
   matrix(
        matern( temp/theta , smoothness = smoothness, scale = scale), 
       nrow = nrow( x1), ncol = nrow( x2) )
}

