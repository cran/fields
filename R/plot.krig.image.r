"plot.krig.image" <-
function (out, main = NA, digits = 4, which = rep(T, 4), graphics.reset = T, 
    ...) 
{
old.par <- par(no.readonly=T)

    if (graphics.reset) {
        on.exit(par(old.par))
        par(xpd = T)
    }
    set.panel(2, 2, T)
    lims <- range(out$fitted.values, out$y)
    if (which[1]) {
        plot(out$fitted.values, out$y, xlim = lims, ylim = lims, 
            ylab = "Observed Values", xlab = "Predicted Values", 
            bty = "n", ...)
        abline(0, 1)
    }
    if (which[2]) {
        image.plot(out$surface, offset = 0.1)
        points(out$x, pch = ".")
    }
    if (which[3]) {
        plot(out$fitted.values, out$residuals, xlim = lims, ylab = "Residuals", 
            xlab = "Predicted\nValues", bty = "n", ...)
        yline(0)
    }
    if (which[4]) {
        look <- list(x = out$grid$x, y = out$grid$y, z = matrix(NA, 
            out$m, out$n))
        look$z[out$indexM] <- out$yM - predict(out, out$xM)
        image.plot(look, xlab = "x", ylab = "y", graphics.reset = F, 
            offset = 0.1)
        contour(out$surface, labex = 0, add = T)
        title("Residuals on surface ", cex = 0.8)
    }
    if (!is.na(main)) 
   mtext(main, cex = 1, outer = T, line = -2)
  #####else      mtext(deparse(out$call), cex = 1, outer = T, line = -2,adj=0)
    invisible()
}
