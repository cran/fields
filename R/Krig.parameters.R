"Krig.parameters" <-
function(obj, mle.calc=TRUE){
# takes fitted estimate and lambda and returns various estiamtes of sigma
# and rho

shat.GCV <- sqrt(sum(obj$weights * 
obj$residuals^2.)/(length(obj$y) - obj$trace))

if(mle.calc) {  rhohat <- sum(out$c * out$yM)/(out$N - out$nt)
shat.MLE <- sqrt(rhohat * out$lambda)}
else{ rhohat<- shat.MLE<- NA}
list( shat.GCV=shat.GCV, shat.MLE= shat.MLE, rhohat= rhohat)
}
