fields.derivative.poly <-
function(x, m,dcoef){

# dimension of x locations
# goal is find partial derivative matrix 

  d<- ncol( x)

  out<- fields.mkpoly( rbind( x[1,]), m)
  ptab<- attr( out, "ptab")

  if( nrow( ptab) != length( dcoef)) { 
       stop(" rows of ptab not equal to length of dcoef")}


  hold<- matrix( NA, ncol=d, nrow= nrow(x))
 
  for( k in 1:d){
     nonzero<- ptab[,k] != 0
     ptemp<- matrix( ptab[nonzero,], ncol=d) 
     
     dtemp<- dcoef[nonzero]
     dtemp<- dtemp *ptemp[,k]
     ptemp[,k]<-  ptemp[,k] -1

     hold[,k]<- fields.evlpoly2(x,dtemp, ptemp)
}

  return(hold)
}

