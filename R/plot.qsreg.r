"plot.qsreg" <-
function (out, pch = "*", main = NA) 
{
    old.par <- par("mfrow", "oma")
    on.exit(par(old.par))
    par(mfrow = c(2, 1))
    plot(out$x, out$y, xlab = "X", ylab = "y", pch=pch)
    orderx <- order(out$x)
    matlines(out$x[orderx], out$fitted.values[orderx, ], lty = 1, 
        col = 1)
    if (nrow(out$cv.grid) > 1) {
        ind <- out$cv.grid[, 3] < 1e+19
        out$cv.grid <- out$cv.grid[ind, ]
        plot(out$cv.grid[, 2], (out$cv.grid[, 3]), xlab = "Effective number of parameters", 
            ylab = "GCV Absolute prediction error ", pch = ".", 
            log = "y")
        title(" Estimated Average Squared Prediction Error", 
            cex = 0.5)
    }
    if (is.na(main)) 
        mtext(deparse(out$call), cex = 1.3, outer = T, line = -2)
    else mtext(main, cex = 1.3, outer = T, line = -2)
}
