  multWendlandGrid <- function( grid.list,center, delta, coef){
     mx<- length( grid.list$x)
     my<- length( grid.list$y)
     dx<- (grid.list$x[mx] - grid.list$x[1])/(mx-1)
     dy<- (grid.list$y[my] - grid.list$y[1])/(my-1)
     centerScaled<- cbind( (center[,1] - grid.list$x[1])/dx + 1,
                       (center[,2] - grid.list$y[1])/dy + 1 )
     deltaX<- delta/dx
     deltaY<- delta/dy
     
     nc<- nrow( center)
     out<-.Fortran( "multWendlandG", PACKAGE="fields",
                 mx=as.integer(mx),
                 my=as.integer(my),
                 deltaX= as.double( deltaX),
                 deltaY= as.double( deltaY),                  
                 nc= as.integer(nc),
                 center=as.double(centerScaled),
                 coef=as.double(coef),
                 h= as.double(matrix(0,mx,my)),
                 flag=as.integer(1)
                )
     if( out$flag!= 0){
       stop("error in multWendlandG FORTRAN")}
     return( out$h)
   }

 predictSurface.fastTps<- function(  object, grid.list, ynew = NULL, 
             Z = NULL, drop.Z = FALSE, just.fixed = FALSE) {
   
    if( ncol(object$x) !=2 | object$args$k != 2){
      stop(" surface function only written for Wendland dim=2, k=2")}
    if( is.list( Z)){
      stop("Z must be passed as a matrix")}
    if (is.null(Z)) {
       drop.Z<- TRUE
    }
    if (!is.null(ynew)) {
        coef.hold <- mKrig.coef(object, ynew)
        c.coef <- coef.hold$c
        d.coef <- coef.hold$d
    }
    else {
        c.coef <- object$c
        d.coef <- object$d
    }
    # fixed part of the model this a polynomial of degree m-1
    # Tmatrix <- fields.mkpoly(xnew, m=object$m)
    #
        xnew<- make.surface.grid( grid.list)
        if (drop.Z | object$nZ == 0) {
            # just evaluate polynomial and not the Z covariate
            temp <- fields.mkpoly(xnew, m = object$m) %*% d.coef[object$ind.drift, ]
        }
        else {
            temp <- cbind(fields.mkpoly(xnew, m = object$m), Z) %*% d.coef
        }
    if (!just.fixed) {
    # add nonparametric part. Covariance basis functions
    # times coefficients.
    temp2<- multWendlandGrid(grid.list, object$knots, delta=object$args$theta, c.coef)
    temp <- temp + temp2
    }
    # add two parts together and coerce to vector
    return( as.surface( grid.list,( temp ) ) )
}

"sim.fastTps.approx" <- function(fastTpsObject, predictionGridlist = NULL,
                            simulationGridlist=NULL, gridRefinement=5, gridExpansion=1,
    M = 1, nx = 40, ny = 40, delta=NULL,  verbose = FALSE ) {
    # create grid if not passed
    if (is.null(predictionGridlist)) {
        predictionGridlist <- fields.x.to.grid(fastTpsObject$x, nx = nx, ny = ny)
        names( predictionGridlist)<- c( "x", "y")
    }
    nx<- length((predictionGridlist$x))
    ny<- length((predictionGridlist$y))
    if ( is.null(simulationGridlist) ) {
        fakePredictionGrid<- cbind( range(predictionGridlist$x), range(predictionGridlist$y))
        simulationGridlist<- makeSimulationGrid( fastTpsObject,fakePredictionGrid ,
                               nx, ny, gridRefinement, gridExpansion)
    }
    nxSimulation<- length(simulationGridlist$x)    
    nySimulation<- length(simulationGridlist$y)
    sigma <- fastTpsObject$sigma.MLE
    rho <- fastTpsObject$rho.MLE
    #
    # set up various sizes of arrays
    nObs <- nrow(fastTpsObject$x)
    if (verbose) {
        cat("nObs, sigma, rho", nObs, sigma, rho, fill = TRUE)
    }
    #
    # set up object for simulating on a grid
    #
#    print( system.time(
    covarianceObject <- wendland.image.cov( 
                            setup = TRUE, grid =simulationGridlist,
                            cov.args=fastTpsObject$args )
#    ))                   
     if (verbose) {
        cat( "dim of full circulant matrix ", dim(covarianceObject$wght), fill = TRUE)
    }
    # output array
    out <- matrix(NA, nx*ny, M ) 
    #
    # find conditional mean field from initial fit
    # don't multiply by sd or add mean if this is
    # a correlation model fit.
    # (these are added at the predict step).
    # from now on all predicted values are on the grid
    # represented by a matrix
   #### hHat <- predict(fastTpsObject, x=predictionGrid)
    hHat<- c(predictSurface.fastTps(fastTpsObject, predictionGridlist)$z)
    # empty image object to hold simulated fields
    hTrue<- c( simulationGridlist, list( z= matrix(NA, nxSimulation,nySimulation)))
    ##########################################################################################
    ### begin the big loop
    ##########################################################################################
    xData<-  fastTpsObject$x
    weightsError<- fastTpsObject$weights
    for (k in 1:M) {
        # simulate full field
        if( verbose){
          cat( k, " ")}
        hTrue$z <- sqrt(rho) * sim.rf(covarianceObject)
        #
        # NOTE: fixed part of model (null space) need not be simulated
        # because the  estimator is unbiased for this part.
        # the variability is still captured because the fixed part
        # is still estimated as part of the predict step below
        #
        #   bilinear interpolation to approximate values at data locations
        #
        hData <- interp.surface(hTrue,xData)
        hPredictionGrid<- c(interp.surface.grid(hTrue, predictionGridlist)$z)
        ySynthetic <- hData + sigma * 1/sqrt(weightsError)* rnorm(nObs)
        # predict at grid using these data
        # and subtract from synthetic 'true' value

        spatialError <- c(predictSurface.fastTps(fastTpsObject, predictionGridlist,
                                               ynew = ySynthetic)$z) - hPredictionGrid
#        spatialError <- predict(fastTpsObject, predictionGrid,
#                                              ynew = ySynthetic) - hPredictionGrid
        # add the error to the actual estimate  (conditional mean)
        out[ , k] <- hHat + spatialError
    }
    return( list( predictionGridlist=predictionGridlist, Ensemble=out, call=match.call()) )
}

makeSimulationGrid<-function( mKrigObject, predictionGrid,
                               nx, ny, gridRefinement, gridExpansion){
        nxSimulation<- nx*gridRefinement*gridExpansion
        nySimulation<- ny*gridRefinement*gridExpansion
        xRange<- range(c(mKrigObject$x[,1], predictionGrid[,1]) )
        yRange<- range(c(mKrigObject$x[,2], predictionGrid[,2]) )
        midpointX<-               (xRange[2] + xRange[1])/2
        midpointY<-               (yRange[2] + yRange[1])/2
        deltaX<-    gridExpansion*(xRange[2] - xRange[1])/2 
        deltaY<-    gridExpansion*(yRange[2] - yRange[1])/2       
        return(
           list( x= seq( midpointX - deltaX, midpointX + deltaX,, nxSimulation),
                 y= seq( midpointY - deltaY, midpointY + deltaY,, nySimulation) )
                )
}

  
