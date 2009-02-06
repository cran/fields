 RadialBasis<- function(d,M,dimension){
# compute the exponent for a thin-plate spline
# based on smoothness and dimension

 p<- 2*M-dimension

 if( p<=0) {
   stop( "M too small for thin plates spline 2m-d >0")}

 if( dimension%%2==0){
   ifelse( d> 1e-14, radbas.constant(M,dimension)*(d^p)* log(d),0)}
 else{
   radbas.constant(M,dimension)*(d^p)}

}

