"Krig" <-
function (x, Y, cov.function = "stationary.cov", lambda = NA, df = NA, 
    cost = 1, knots=NA, weights = NULL, m = 2, 
    nstep.cv = 80, scale.type = "user", x.center = rep(0, ncol(x)), 
    x.scale = rep(1, ncol(x)), rho = NA, sigma2 = NA, method = "GCV", 
     verbose = FALSE, mean.obj = NA, 
    sd.obj = NA, 
    null.function = fields.mkpoly, 
    offset = 0, outputcall = NULL, cov.args = NULL, na.rm = TRUE,
    chol.args=NULL, give.warnings=TRUE,
    ...) 
# NOTES 
# the verbose switch prints many intermediate steps as an aid in debugging.
#
{

#
# create output list 
    out <- list()

###########################################################
#  First series of steps simply store pieces of the passed 
#    information to output list (i.e. the Krig object)
##########################################################

  if( is.null( outputcall)){
       out$call<- match.call()}
    else{
        out$call <- outputcall}
#
# save covariance function as it name 
#
   if (!is.character(cov.function)) {
        if (is.function(cov.function)) {
            out$cov.function.name <- as.character(substitute(cov.function))}
        
    }
    else{
      out$cov.function.name <- cov.function }

#
# logical to indicate if the "C" argument is present in this function
    C.arg.missing<- all( names( formals( get( out$cov.function.name)))!="C")
    if( C.arg.missing) stop("Need to have C argument in covariance function 
                                 see exp.cov.simple as an example")

#
# save parameters values possibly passed to the covariance function
    if (!is.null(cov.args)) 
        out$args <- cov.args
    else out$args <- list(...)
#
# default values for Cholesky decomposition, these are important
# for sparse matrix decompositions used in Krig.engine.fixed. 

    if( is.null( chol.args)) {
        out$chol.args<- list( pivot= FALSE)}
    else{
        out$chol.args<- chol.args}


#
#
# save function that creates null space.
# for fields.mkpoly  the order of polynomial will be (m-1)

    if( m < 1) { 
       stop("m needs to be 1 or greater" )}
    out$make.tmatrix <- null.function
    out$m <- m 
#
# the offset is the effective number of parameters used in the GCV 
# calculations
    out$offset <- offset

#
# the cost is the multiplier applied to the GCV eff.df
# sigma2 is error variance and rho the multiplier for covariance
    out$cost <- cost
    out$sigma2<- sigma2
    out$rho<- rho
#
# correlation model information
#
    out$mean.obj<- mean.obj
    out$sd.obj<- sd.obj
    out$correlation.model <- !(is.na(mean.obj[1])&is.na( sd.obj[1]))

#
# transformation info
   out$scale.type<- scale.type
   out$x.center<- x.center
   out$x.scale<- x.scale

#
# verbose block
    if (verbose) {
        cat(" Local cov function arguments ", fill = TRUE)
        print(out$args)
        cat(" covariance function used is : ", fill = TRUE)
        print(out$cov.function.name)
    }

###############################################################
# Begin modifications and transformations of input information
###############################################################
 
# various checks on x and  Y including removal of NAs in Y
   out2<- Krig.check.xY( x,Y, weights, na.rm, verbose=verbose)
   out<- c( out, out2)


# transform to correlation model (if appropriate)
# find replicates and collapse to means and pool variances.
# Transform unique x locations and knots. 

 
   if( out$correlation.model){
               out$y<- Krig.cor.Y(out, verbose=verbose)}

   out2<- Krig.transform.xY(out,knots, verbose=FALSE)
   out<- c( out, out2)


# NOTE: knots have been transformed after this step


#############################################################
#  Figure out what to do 
#############################################################

#
# determine the method for finding lambda 
#  Note order 

    out$method<- method

    if (!is.na(lambda)  ){
# this indicates lambda has been supplied and leads to 
# the cholesky type computational approaches. 
        out$method <- "user"
        out$lambda<- lambda
    }

    if (!is.na(rho) & !is.na(sigma2)) {
        out$method <- "user"
        out$lambda <- sigma2/rho
    }

#
# NOTE: method="user" means that a value of lambda has been supplied
#        and so GCV etc to determine lambda is not needed. 
#     
   out$fixed.model<- (out$method=="user")

# set lambda.est matrix to NA because no estimates are found
# (see alternative in gcv block)

    if( out$fixed.model) { out$lambda.est<- NA}

#
# verbose block 
   if (verbose){ 
        cat("lambda, fixed? ", lambda, out$fixed.model, fill = TRUE)}



########################################################
#   Do the intensive linear algebra to find the solution
########################################################

# this is where all the heavy lifting happens and are termed the 
# engines. 
#
# Note that all the information is passed as a list
# incl,uding arguments to the cholesky decomposition 
# used within Krig.engine.fixed
#
# The results are saved in the component _matrices_ 
#
# if method=="user" then just evaluate at single lambda
#  _fixed_ here means a fixed lambda
#
# For fixed lambda the decompostions with and without knots
# are surprisingly similar and so are in one function
#

  if( out$fixed.model){
       out$matrices<-  Krig.engine.fixed( out, verbose=verbose)
   # can't find the trace of A matrix in fixed lambda case. 
       out$eff.df<- NA
  }
#
# alternative are 
# matrix decompositions suitable for 
# evaluation at many lambdas to facilitate GCV/REML estimates  etc. 
#
  if( !out$fixed.model){
    if( out$knot.model){
           out$matrices <- Krig.engine.knots( out, verbose=verbose)
           out$pure.ss <- out$matrices$pure.ss}

    else{
           out$matrices<- Krig.engine.default( out, verbose=verbose)
    }
  } 
#
# store basic information about decompositions

 out$nt<- out$matrices$nt # dim of null space
 out$np<- out$matrices$np # number of basis functions
 out$decomp<- out$matrices$decomp # type of decomposition see Krig.coef


#################################################
# Do GCV and REML search over lambda if not fixed
#################################################
  if( !out$fixed.model){

    if(verbose){ cat("call to gcv.Krig", fill=TRUE)}

      gcv.out <- gcv.Krig(out, nstep.cv = nstep.cv, verbose = verbose, 
            cost = out$cost, offset = out$offset, give.warnings=give.warnings)
      out$gcv.grid <- gcv.out$gcv.grid
#
#  a handy summary table of the search results
      out$lambda.est <- gcv.out$lambda.est
#
# verbose block
      if (verbose) {
           cat("returned GCV and REML grid search", fill=TRUE)
           print(out$gcv.grid)
       }
#
# assign the preferred lambda either from GCV/REML/MSE or the user value
#
       out$lambda <- gcv.out$lambda.est[out$method, 1]
       out$eff.df<- out$lambda.est[out$method, 2] 
        
       if (verbose) {
            cat("trace of A", fill = TRUE)
            print(out$eff.df)
        }
    }
# end GCV/REML block 

##########################################
# find coefficients at prefered lambda 
# and evaluate the solution at observations
##########################################
#   pass replicate group means -- noneed to recalculate these. 

    out2 <- Krig.coef(out, yM= out$yM)
    out<- c( out, out2)
   
#
# fitted values and residuals and predicted values on null space (fixed 
# effects). But be sure to do this at the nonmissing x's
#
    out$fitted.values <- predict.Krig(out, out$x, 
                            eval.correlation.model = FALSE)
    out$residuals <- out$y - out$fitted.values
    out$fitted.values.null <- as.matrix(
                out$make.tmatrix(out$x, m)) %*% out$d 
#
# verbose block
    if (verbose) {
        cat("residuals", out$residuals, fill = TRUE)
    }
#

# find various estimates of sigma and rho 

      out2<-Krig.parameters(out)
      out<- c( out, out2)
#
# assign the "best" model as a default choice 
# either use the user supplied values or the results from 
#  optimization
#
      if(out$method=="user"){
           out$best.model <- c(out$lambda, out$sigma2, out$rho)}
      else{
          out$best.model <- c(out$lambda, out$shat.MLE**2, out$rhohat)}

# Note: values in best.model are used in subsquent functions as the choice 
# for these parameters!


# set class 
    class(out) <- c("Krig")

    return(out)
}

