"plot.qsreg" <-
function (out, pch = "*", main = NA) 
{
    old.par <- par("mfrow", "oma")
    on.exit(par(old.par))
    set.panel( 2, 2, relax=TRUE)
    plot(out$x, out$y, xlab = "X", ylab = "y", pch = pch)
    orderx <- order(out$x)
temp<-  out$fitted.values[,c(out$ind.cv, out$ind.cv.ps)]
matlines(out$x[orderx], temp[orderx,], lty = 1,  col = c(1,2))
##
# residual plot
#
matplot(out$x, 
qsreg.psi( out$residuals[,c(out$ind.cv, out$ind.cv.ps)],  out$alpha, out$sc)
, col=c(1,2), pch="o", ylab="Pseudo residuals", xlab="X")

yline(0)
   

    if (nrow(out$cv.grid) > 1) {
        ind <- out$cv.grid[, 3] < 1e+19
        out$cv.grid <- out$cv.grid[ind, ]
        matplot(out$cv.grid[, 2], cbind(out$cv.grid[, 3], 
out$cv.grid[,6]), xlab = "Effective number of parameters", 
            ylab = "Log CV Rho function ", 
            log = "y", type="l", col=c(1,2))
xline( out$cv.grid[out$ind.cv, 2], col=1)
xline( out$cv.grid[out$ind.cv.ps, 2], col=2)

        title(" CV curves", 
            cex = 0.5)
    }

bplot(  
qsreg.psi( out$residuals[,c(out$ind.cv, out$ind.cv.ps)],  out$alpha, out$sc), 
labels=c("CV", "CV pseudo")
)

yline( 0, col=2)

    if (is.na(main)) 
        mtext(deparse(out$call), cex = 1.3, outer = TRUE, line = -2)
    else mtext(main, cex = 1.3, outer = TRUE, line = -2)
}
