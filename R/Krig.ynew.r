# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html

"Krig.ynew" <-
function(out, y=NULL, yM=NULL)
{
	#
	# calculates the collapsed y (weighted) mean vector based on the
	# X matrix and weights from the out object. 
	# or just passes through the collapsed mean data if passed. 
        #
	#
	# If there are no replicated obs. then return the full vector
	# pure error ss is zero 
	#

		shat.rep <- NA
		shat.pure.error <- NA

		pure.ss <- 0

# if no y's are given then it is assumed that one should use the
# yM from the original data used to create the Krig object

        if( is.null(yM)&is.null( y) ) {
                        yM<- out$yM}


#
# case when yM is passed no calculations are needed
#
        if( !is.null(yM) ){
                return(list(yM = as.matrix(yM),
                            shat.rep = NA, shat.pure.error = NA,
                        pure.ss = 0))
        }


#
# no reps case
#
	if(length(unique(out$rep.info)) == out$N) {
		return(list(yM = as.matrix(y),
                            shat.rep = NA, shat.pure.error = NA,
			pure.ss = 0))
	}
# 
#  check that y is the right length
#

	if(length(y) != out$N) {
		stop(" the new y vector is the wrong length!")
	}
#
# case when full y data is passed and replicate means need to be found
#
	if(length(unique(out$rep.info)) < out$N) {
           
#
# calculate means by pooling Replicated obseravations but use the
# the right weighting. 
#
		rep.info.aov <- fast.1way(out$rep.info, y, out$weights)[
			c("means", "MSE", "SSE")]
		shat.pure.error <- sqrt(rep.info.aov$MSE)
		shat.rep <- shat.pure.error
		return(list(yM = rep.info.aov$means, shat.rep = shat.rep, 
			shat.pure.error = shat.pure.error, pure.ss = 
			rep.info.aov$SSE))
	}


}
