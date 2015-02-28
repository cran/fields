printGCVWarnings<- function( Table, method="all"){
 ind<- Table$Warning
 if( method == "all"){
 	kIndex<- 1:6
 }
 else{
    kIndex<- match( method,c("GCV",
                          "GCV.model",
                          "GCV.one",
                          "RMSE",
                          "pure error",
                         "REML")
                   )
      }
 methodList<- c(
 "(GCV) Generalized Cross-Validation ",
 "(GCV.model) Generalized Cross-Validation on replicate means ",
 "(GCV.one) Generalized Cross-Validation on individual observations ",
 "(RMSE) Matching estimate of sigma to supplied rmse ",
 "Matching estimate of sigma to that from replicated observations",
 "(REML) Restricted maximum likelihood "
  )
 if( any( ind[kIndex])){
 	cat("Warning: ", fill=TRUE)
 	cat("Grid searches over lambda (nugget and sill variances) with  minima at the endpoints: ", fill=TRUE) }
 for( k in kIndex){
 if( ind[k]){
    whichEnd<- ifelse(Table[k,2],"left","right")
 	cat( " ", methodList[k], fill =TRUE)
 	cat( "   minimum at ", whichEnd, "endpoint",
 	                " lambda  = ", Table[k,6] ,
    "(eff. df=", Table[k,7] , ")", fill = TRUE )
 	     }
 }	     
 
}