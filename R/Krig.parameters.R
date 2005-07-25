"Krig.parameters" <-
function (obj, mle.calc = obj$mle.calc) 
{
    shat.GCV <- sqrt(sum(obj$weights * obj$residuals^2)/(length(obj$y) - 
        obj$eff.df))
    if (mle.calc) {
        rhohat <- sum(obj$c * obj$yM)/(obj$N - obj$nt)
        shat.MLE <- sqrt(rhohat * obj$lambda)
    }
    else {
        rhohat <- shat.MLE <- NA
    }
    list(shat.GCV = shat.GCV, shat.MLE = shat.MLE, rhohat = rhohat)
}

