"gcv.Krig" <-
function (out, lambda.grid = NA, cost = 1, nstep.cv = 80, rmse = NA, 
    verbose = FALSE, tol = 1e-05, offset = 0, y = NULL, 
give.warnings=TRUE) 
{
# set up some temporary variables. 
    nt <- out$nt
    np <- out$np
    N <- out$N

# D are the eigenvalues of the "covariance" matrix form which most 
# gnarly things are compuated. 

    D <- out$matrices$D

# are obs supplied? if not use object 
#
    if (is.null(y)) {
        u <- out$matrices$u
        shat.pure.error <- out$shat.pure.error
        pure.ss <- out$pure.ss
    }
    else {
        out2 <- Krig.updateY(out, y)
        u <- out2$u
        shat.pure.error <- out2$shat.pure.error
        pure.ss <- out2$pure.ss
    }
#
# figure out a good range for lambda
# the most rational is to base it on the model effective degrees of freedom
#
    if (is.na(lambda.grid[1])) {
#
# create some pretty values for the lambda grid
#
# set up grid initially in terms of effective degrees of freedom 
        temp.df<- seq(nt, (np - offset) * 0.95,,nstep.cv)
#
# offset the first value because this is infinity!
# 
        temp.df[1] <- temp.df[1] + .001
        for( k in 1: nstep.cv){
                   lambda.grid[k] <- Krig.df.to.lambda(temp.df[k], D)
        }
    }

#
# sort lambda.grid to be consistent with old versions of Krig
#

lambda.grid<- sort( lambda.grid) 

# set up some variables for the coarse grid search

    nl <- length(lambda.grid)
    nd <- length(D)
    V<- V.model<- V.one<-lplike<-trA<-shat<- rep(NA, nl)
    Dl <- rep(NA, nd)
 
#
# this is a small "Krig" object to make it easier to pass these pieces to 
# different gcv functions. See f.extra in golden.section search.
#
    info <- list(matrices = list(D = D, u = u), N = N, nt = nt, 
        cost = cost, pure.ss = pure.ss, shat.pure.error = shat.pure.error, 
        offset = offset)
#
# this loop is a bit inefficient because many quantities are being 
# recalculated but it makes it clear from where these different
# versions are coming. 
#
   for (k in 1:nl) {
        V[k] <- Krig.fgcv(lambda.grid[k], info)
        V.one[k] <- Krig.fgcv.one(lambda.grid[k], info)
        V.model[k]<- Krig.fgcv.model(lambda.grid[k], info)
        lplike[k]<- Krig.flplike(lambda.grid[k], info)
        shat[k]<-  sqrt(Krig.fs2hat(lambda.grid[k], info))
        trA[k] <- Krig.ftrace(lambda.grid[k],D)
    }
 
# combine these values into a data frame. This is the information 
# used in the 3rd plot in the plot.Krig figure. 

    gcv.grid <- cbind(lambda.grid, trA, V, V.one, V.model, shat, lplike)
    gcv.grid <- as.data.frame(gcv.grid)
    names(gcv.grid) <- c("lambda", "trA", "GCV", "GCV.one", "GCV.model", 
        "shat","-Log Profile")
    if (verbose) 
        print(info)
    if (verbose) 
        print(gcv.grid)
#
# at this point one could simply interpolate the coarse search 
#( using sreg of course!) over the lambda
# grid and find the minimum. A more deliberate way to do is a refined 
# optimization using a godlen section search. 
#

# create the output summary matrix

lambda.est <- matrix(ncol =4 , nrow = 6, 
dimnames = list(
  c("GCV", "GCV.model", "GCV.one", "RMSE", "pure error", "-Log Profile" ), 
   c("lambda", "trA", "GCV", "shat"))
)
#
# Now start filling in lamnda.est with the refined estimates
# using various versions of GCV and other approaches. 
#
# refine the traditional GCV estimate
# If there are replicates it is leave-one-replicate-group- out. 
#
    lambda.est[1, 1] <- Krig.find.gcvmin(info, lambda.grid, gcv.grid$GCV, 
        Krig.fgcv, tol = tol, verbose = verbose)
#
# Next find lambda so that the estimate of sigma**2 matches the 
# the estimate based on the replicated points. 

   if (!is.na(shat.pure.error)) {
        temp <- gcv.grid$GCV.model
        lambda.est[2, 1] <- Krig.find.gcvmin(info, lambda.grid, 
            temp, Krig.fgcv.model, tol = tol, verbose = verbose)
    }
#
# GCV using really leave-one-out. 
#
    lambda.est[3, 1] <- Krig.find.gcvmin(info, lambda.grid, gcv.grid$GCV.one, 
        Krig.fgcv.one, tol = tol, verbose = verbose)

# Profile likelihood
lambda.est[6, 1] <- Krig.find.gcvmin(info, lambda.grid, 
 gcv.grid$"-Log Profile",
        Krig.flplike, tol = tol, verbose = verbose)
if( verbose){
cat( " mle estimate", lambda.est[6, 1], fill=TRUE)}



#
# If rmse or if shat.pure.error can is available find the lambda that 
# yeilds an esimate that mathces this value. 
# This is a zero-finding task so we use a simple algorithm that 
# determines zero crossing. 
#
    lambda.rmse <- NA
    lambda.pure.error <- NA
# 
# work on matching the passed rmse value
    if (!is.na(rmse)) {
        if (all(gcv.grid$shat < rmse) | all(gcv.grid$shat > rmse)) {
            guess <- NA
        }
        else {
            guess <- max(gcv.grid$lambda[gcv.grid$shat < rmse])
        }
        if (verbose) {
            print(rmse)
            print(guess)
        }
        if (!is.na(guess)) {
            lambda.rmse <- find.upcross(Krig.fs2hat, info, upcross.level = rmse^2, 
                guess = guess, tol = tol * rmse^2)
            lambda.est[4, 1] <- lambda.rmse
        }
        else {
            if( give.warnings){
             warning("Value of rmse is outside possible range")}
        }
    }
#
# work on matching shat.pure.error

   if (!is.na(shat.pure.error)) {
        if (all(gcv.grid$shat < shat.pure.error) | all(gcv.grid$shat > 
            shat.pure.error)) {
            guess <- NA
        }
        else {
            guess <- max(gcv.grid$lambda[gcv.grid$shat < shat.pure.error])
        }
        if (!is.na(guess) & (guess != -Inf)) {
            lambda.pure.error <- find.upcross(Krig.fs2hat, info, 
                upcross.level = shat.pure.error^2, guess = guess, 
                tol = tol * shat.pure.error^2)
            lambda.est[5, 1] <- lambda.pure.error
        }
        else {
            if( give.warnings){
             warning("Value of pure error estimate  is outside possible range")}
        }
    }
if( verbose) {
cat( " Done with refining", fill=TRUE)}

#
# fill in the rest of the lambda.est summary matrix
# This is reporting the values of various criterion at the 
# lambda found by GCV traditional ( leave-one/group-out)
#
    for (k in 1:5) {
        lam <- lambda.est[k, 1]
        if (!is.na(lam)) {
            lambda.est[k, 2] <- Krig.ftrace(lam, D)
            if (k == 1 | k > 3) {
                lambda.est[k, 3] <- Krig.fgcv(lam, info)
            }
            if (k == 2) {
                lambda.est[k, 3] <- Krig.fgcv.model(lam, info)
            }
            if (k == 3) {
                lambda.est[k, 3] <- Krig.fgcv.one(lam, info)
            }
            
            lambda.est[k, 4] <- sqrt(Krig.fs2hat(lam, info))
        }
    }
#
# fill in some information about ML estimate of lambda

      lam.ml<-  lambda.est[6,1]
      lambda.est[6,2] <-  Krig.ftrace(lam.ml, D)
      lambda.est[6,3] <-  Krig.fgcv(lam.ml, info)
      lambda.est[6,4] <- sqrt(Krig.fs2hat(lam.ml, info))
#
# and now the output list
#
    list(gcv.grid = gcv.grid, lambda.est = lambda.est, 
            lambda.best = lambda.est[1,1])
}
