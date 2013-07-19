# fields, Tools for spatial data
# Copyright 2004-2011, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
"sim.mKrig.approx" <- function(mKrigObject, predictionGrid = NULL,
                            simulationGridlist=NULL, gridRefinement=5, gridExpansion=1,
    M = 1, nx = 40, ny = 40,  verbose = FALSE ) {
    # create grid if not passed
    if (is.null(predictionGrid)) {
        predictionGridlist <- fields.x.to.grid(mKrigObject$x, nx = nx, ny = ny)
        predictionGrid<- make.surface.grid( predictionGridlist)
    }
    if (is.null(simulationGridlist)) {
        nxSimulation<- nx*gridRefinement*gridExpansion
        nySimulation<- ny*gridRefinement*gridExpansion
        xRange<- range(c(mKrigObject$x[,1], predictionGrid[,1]) )
        yRange<- range(c(mKrigObject$x[,2], predictionGrid[,2]) )
        midpointX<-               (xRange[2] + xRange[1])/2
        midpointY<-               (yRange[2] + yRange[1])/2
        deltaX<-    gridExpansion*(xRange[2] - xRange[1])/2 
        deltaY<-    gridExpansion*(yRange[2] - yRange[1])/2       
        simulationGridlist <- list( x= seq( midpointX - deltaX,  midpointX + deltaX,, nxSimulation),
                                    y= seq( midpointY - deltaY,  midpointY + deltaY,, nySimulation) )
    }
     if (verbose) {
        cat("predictionGrid summary stats",  fill = TRUE)
        print( stats(predictionGrid))
        print(stats(simulationGridlist), fill = TRUE)
        
    }
        sigma <- mKrigObject$sigma.MLE
        rho <- mKrigObject$rho.MLE
    #
    # set up various sizes of arrays
    nObs <- nrow(mKrigObject$x)
    if (verbose) {
        cat("nObs, sigma, rho", nObs, sigma, rho, fill = TRUE)
    }
    #
    # set up object for simulating on a grid
    #
    covarianceObject <- stationary.image.cov( 
                            setup = TRUE, grid =simulationGridlist,
                            cov.function=mKrigObject$cov.function,  mKrigObject$args )
     if (verbose) {
        cat( "dim of full circulant matrix ", dim(covarianceObject$wght), fill = TRUE)
    }
    # output array
    out <- matrix(NA, nrow(predictionGrid), M ) 
    #
    # find conditional mean field from initial fit
    # don't multiply by sd or add mean if this is
    # a correlation model fit.
    # (these are added at the predict step).
    # from now on all predicted values are on the grid
    # represented by a matrix
    hHat <- predict(mKrigObject, x=predictionGrid)
    # empty image object to hold simulated fields
    hTrue<- c( simulationGridlist, list( z= matrix(NA, nxSimulation,nySimulation)))
    ##########################################################################################
    ### begin the big loop
    ##########################################################################################
    xData<-  mKrigObject$x
    weightsError<- mKrigObject$weights
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
        hPredictionGrid<- interp.surface(hTrue, predictionGrid)
        ySynthetic <- hData + sigma * 1/sqrt(weightsError)* rnorm(nObs)
        # predict at grid using these data
        # and subtract from synthetic 'true' value
        spatialError <- predict(mKrigObject, xnew=predictionGrid, y = ySynthetic) - hPredictionGrid
        # add the error to the actual estimate  (conditional mean)
        out[ , k] <- hHat + spatialError
    }
    return( list( predictionGrid=predictionGrid, Ensemble=out, call=match.call()) )
}
