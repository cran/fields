"plot.Krig"<-
function(out, main = NA, digits = 4, which = c(T, T, T, T), graphics.reset = T, 
	...)
{
	old.par <- par("mfrow", "oma")
	if(graphics.reset) {
		on.exit(par(old.par))
		par(xpd = T)
	}
	set.panel(2, 2, T)
	fitted.values <- predict(out)
	std.residuals <- (out$residuals * sqrt(out$weights))/out$shat.GCV
	if(which[1]) {
		temp <- summary(out)
		plot(fitted.values, out$y, ylab = "Y", xlab = 
			" predicted values", bty = "n", ...)
		abline(0, 1)
		hold <- par("usr")
		text(hold[1], hold[4], paste(" R**2 = ", format(round(100 * 
			temp$covariance, 2)), "%", sep = ""), cex = 0.8, adj = 
			0)
	}
	if(which[2]) {
		plot(fitted.values, std.residuals, ylab = "(STD) residuals", 
			xlab = " predicted values", bty = "n", ...)
		yline(0)
		hold <- par("usr")
		text(hold[1], hold[4], paste(" RMSE =", format(signif(sqrt(sum(
			out$residuals^2)/(temp$num.observation - temp$enp)), 
			digits))), cex = 0.8, adj = 0)
	}
	if(which[3]) {
		if(nrow(out$gcv.grid) > 1) {
# trim off + infinity due to pole in the denominator of GCV function
#with cost
			ind <- out$gcv.grid[, 3] < 1e+019
			out$gcv.grid <- out$gcv.grid[ind,  ]
			yr <- range(unlist(out$gcv.grid[, 3:5]))
			plot(out$gcv.grid[, 2], out$gcv.grid[, 3], xlab = 
				"Eff. number of parameters", ylab = 
				" GCV function", bty = "n", ylim = yr, log = 
				"y", ...)
			lines(out$gcv.grid[, 2], out$gcv.grid[, 4], lty = 2)
			lines(out$gcv.grid[, 2], out$gcv.grid[, 5], lty = 1)
			xline(out$eff.df)	#			hold <- par("usr")
#			text(out$eff.df, hold[4], paste(" Eff. df. =", 
#format(
#				round(out$eff.df, 1)), "\n Res. df. =", 
#format(
#				round(temp$num.observation - temp$enp, 
#1))), 
#				cex = 0.8, adj = 0)
			title("GCV-points , solid- GCV model,\ndashed- GCV one",
				cex = 0.6)
		}
	}
	if(which[4]) {
		hist(std.residuals)
	}
	if(is.na(main))
		mtext(deparse(out$call), cex = 1.3, outer = T, line = -2)
	else mtext(main, cex = 1.3, outer = T, line = -2)
}
