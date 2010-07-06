RadialBasis <- function(d, M, dimension, derivative=0) {
    # compute the exponent for a thin-plate spline
    # based on smoothness and dimension
    p <- 2 * M - dimension
    if (p <= 0) {
        stop("M too small for thin plates spline, need: 2m-d >0")
    }
    if( (p-1 <0) & (derivative >0) ) {
      stop("M is too small for derivatives, need: 2m-d < 1")}   
    if( derivative==0){
      if (dimension%%2 == 0) {
        # factor of 2 from the log term
        ifelse(d > 1e-14, radbas.constant(M, dimension) * 
          (d^p) * log(d), 0)}
      else {
        radbas.constant(M, dimension) * (d^p)}}
## find derivative    
    else{

        if (dimension%%2 == 0) {
        # factor of 2 from the log term
          ifelse(d > 1e-14,
                 radbas.constant(M, dimension) * (d^(p-1))* (p* log(d) +1) , 0)}
        else {
         con<-  radbas.constant(M, dimension)*p
         con * (d^(p-1))}}
##### should not get here!    
}
