#
# fields  is a package for analysis of spatial data written for
# the R software environment.
# Copyright (C) 2024 Colorado School of Mines
# 1500 Illinois St., Golden, CO 80401
# Contact: Douglas Nychka,  douglasnychka@gmail.com,
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2
##END HEADER
"summary.sreg" <- function(object, digits = 4, ...) {
    x <- object
    if (length(x$lambda) > 1) {
        stop("Can't do a summary on an object with a grid of smoothing\nparameters")
    }
    summary <- list(call = x$call, num.observation = length(x$residuals), 
        enp = x$trace, nt = x$nt, res.quantile = quantile(x$residuals, 
            seq(0, 1, 0.25)), tauHat.GCV = x$tauHat.GCV, m = x$m, 
        lambda = x$lambda, cost = x$cost, num.uniq = length(x$y), 
        np = x$np, method = x$method, lambda.est = x$lambda.est[!is.na(x$lambda.est[, 
            1]), ], tauHat.pure.error = x$tauHat.pure.error)
    class(summary) <- "summary.sreg"
    summary$covariance <- cor(x$fitted.values * sqrt(x$weights), 
        (x$y) * sqrt(x$weights))^2
    hold <- (sum((x$y - mean(x$y))^2) - sum(x$residuals^2))/(sum((x$y - 
        mean(x$y))^2))
    summary$adjr2 <- 1 - ((length(x$residuals) - 1)/(length(x$residuals) - 
        x$eff.df)) * (1 - hold)
    summary$digits <- digits
    summary$sum.gcv.lambda <- summaryGCV.sreg(x, x$lambda)
    summary
}
