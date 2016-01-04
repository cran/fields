# fields, Tools for spatial data
# Copyright 2004-2007, Institute for Mathematics Applied Geosciences
# University Corporation for Atmospheric Research
# Licensed under the GPL -- www.gpl.org/licenses/gpl.html
supportsArg = function(fun=stationary.cov, arg) {
  
  if(is.null(fun)) {
    #set fun to the default covariance function if not specified
    fun = stationary.cov
  }
  
  argNames = names(as.list(args(fun)))
  return(any(argNames == arg))
}
