\name{fields internal }
\alias{[.spatial.design}
\alias{COR}
\alias{D4transform.image}
\alias{FIELDS.version}
\alias{Krig.coef}
\alias{Krig.df.to.lambda}
\alias{Krig.discretize}
\alias{Krig.fdf}
\alias{Krig.fgcv}
\alias{Krig.fgcv.model}
\alias{Krig.fgcv.one}
\alias{Krig.find.gcvmin}
\alias{Krig.flplike}
\alias{Krig.fs2hat}
\alias{Krig.ftrace}
\alias{Krig.parameters}
\alias{Krig.replicates}
\alias{Krig.updateY}
\alias{Krig.ynew}
\alias{WD4}
\alias{WD42d}
\alias{WD42di}
\alias{WD4i}
\alias{WQS}
\alias{WQS.T}
\alias{WQS.basis}
\alias{WQS2d}
\alias{WQS2di}
\alias{WQSdi}
\alias{WQSi}
\alias{WQSi.T}
\alias{WQS2d.cylinder}
\alias{WQS2di.cylinder}
\alias{WQSi.periodic} 
\alias{WQSi.periodic.T} 
\alias{WQS.periodic} 
\alias{WQS.periodic.basis}
\alias{WQS.periodic.T} 
\alias{bisection.search}
\alias{bplot.obj}
\alias{cat.matrix}
\alias{cat.to.list}
\alias{ceiling2}
\alias{conjugate.gradient}
\alias{describe}
\alias{dyadic.2check}
\alias{dyadic.check}
\alias{describe.bplot}
\alias{double.exp}
\alias{draw.bplot}
\alias{draw.bplot.obj}
\alias{exp.cov.S}
\alias{exp.earth.cov}
\alias{fast.1way}
\alias{find.upcross}
\alias{gauss.cov}
\alias{gcv.sreg}
\alias{golden.section.search}
\alias{grey.level}
\alias{image.plot.info}
\alias{image.plot.plt}
\alias{image.smooth.setup}
\alias{in.poly}
\alias{krig.image.parameters}
\alias{lonlat2xy}
\alias{make.surface.grid}
\alias{make.tmatrix}
\alias{make.tmatrix.F}
\alias{make.tmatrix.S}
\alias{matern}
\alias{minimax.crit}
\alias{periodic.cov.1d}
\alias{periodic.cov.cyl}
\alias{periodic.plane.3d}
\alias{plot.coef}
\alias{plot.krig.image}
\alias{plot.qsreg}
\alias{plot.sim.krig.image}
\alias{plot.spatial.design}
\alias{plot.vgram.matrix}
\alias{predict.interp.surface}
\alias{predict.krig.image}
\alias{predict.qsreg}
\alias{predict.sreg}
\alias{predict.Tps}
\alias{print.krig.image}
\alias{print.qsreg}
\alias{print.spatial.design}
\alias{print.sreg}
\alias{print.summary.Krig}
\alias{print.summary.krig.image}
\alias{print.summary.spatial.design}
\alias{print.summary.sreg}
\alias{qr.q2ty}
\alias{qr.yq2}
\alias{qsreg.fit}
\alias{qsreg.psi}
\alias{qsreg.rho}
\alias{qsreg.trace}
\alias{rad.cov}
\alias{rad.cov.F}
\alias{rad.cov.S}
\alias{rad.covC}
\alias{rad.image.cov}
\alias{radbas.constant}
\alias{rdist.F}
\alias{rdist.S}
\alias{replace.args.function}
\alias{sim.krig.image}
\alias{sim.rf.default}
\alias{sreg.df.to.lambda}
\alias{sreg.fdf}
\alias{sreg.fgcv}
\alias{sreg.fgcv.model}
\alias{sreg.fgcv.one}
\alias{sreg.fit}
\alias{sreg.fs2hat}
\alias{sreg.trace}
\alias{stats.bplot}
\alias{stats.sim.krig.image}
\alias{summary.gcv.Krig}
\alias{summary.gcv.sreg}
\alias{summary.krig.image}
\alias{summary.qsreg}
\alias{summary.spatial.design}
\alias{summary.sreg}
\alias{surface}
\alias{surface.default}
\alias{surface.krig.image}
\alias{surface.surface}
\alias{test.data}
\alias{test.data2}
\alias{unscale}
\alias{variogram.matrix}
\alias{world.dat}
\title{
  Fields internal and secondary functions 
}
\description{
Listed below are supporting fucntions for the major methods in fields. 
}
\usage{

COR (dat) 
 
D4transform.image (x, inv = FALSE, transpose = FALSE, cut.min = 8) 
 
FIELDS.version  

Krig.coef (out, u = NULL, lambda = NULL, yM = NULL) 
Krig.df.to.lambda (df, D, guess = 1, tol = 1e-05) 
Krig.discretize (x, m = 64, n = 64, xg = NULL, yg = NULL, grid 
    = NULL,  expand = c(1, 1)) 
Krig.fdf (llam, info) 
Krig.fgcv (lam, obj) 
Krig.fgcv.model (lam, obj) 
Krig.fgcv.one (lam, obj) 
Krig.find.gcvmin (info, lambda.grid, gcv, gcv.fun, tol, verbose =FALSE, give.warnings = TRUE)
Krig.flplike (lam, obj)
Krig.fs2hat (lam, obj) 
Krig.ftrace (lam, D) 
Krig.parameters (obj, mle.calc=TRUE)
Krig.replicates (out, verbose = FALSE) 
Krig.updateY (out, Y, verbose = FALSE) 
Krig.ynew (out, ynew, verbose = FALSE) 
 

WD4 (x) 
WD42d (x) 
WD42di (x) 
WD4i (x) 
 
WQS (x) 
WQS.T (x) 
WQS.basis (N, cut.n = 8) 
WQS2d (x, transpose = FALSE) 
WQS2di (x, transpose = FALSE) 
WQSdi (x, transpose = FALSE) 
WQSi (x) 
WQSi.T (x) 
WQS.periodic(x)
WQS.periodic.T(x)
WQS.periodic.basis(N, cut.n = 8)
WQS.periodic.T(x)
WQS2d.cylinder(x, transpose = FALSE, byX=T)
WQS2di.cylinder(x, transpose = FALSE, byX=T)
WQSi.periodic(x)
WQSi.periodic.T(x)


bisection.search (x1, x2, f, tol = 1e-07, niter = 25, f.extra = 
        NA, upcross.level = 0) 
 
bplot.obj (data, pos = NA, width = NULL, labels = NULL, las=NULL
            ,add = FALSE, space = 0.25, sort.names = FALSE, xlab = "", ylab = "", 
              label.cex = 1, xaxt = "n", outlier = TRUE, horizontal = FALSE, ...) 
 
cat.matrix (mat, digits = 8) 
 
cat.to.list (x, a) 
 
ceiling2 (m) 
 
conjugate.gradient (b, multAx, start = NULL, tol = 1e-05, kmax = 
             25, verbose = TRUE,    ...) 
 
describe (x) 
describe.bplot (temp, style = "tukey", outlier = TRUE) 

double.exp(x)

dyadic.2check( m,n,cut.p=2)
dyadic.check( n,cut.p=2) 

draw.bplot (temp, width, xpos, outlier = TRUE, style = "tukey") 
draw.bplot.obj (obj, width, xpos, outlier = TRUE, horizontal = FALSE) 
 
exp.cov.S (x1, x2, theta = 1) 
exp.earth.cov (x1, x2, theta = 1) 
 
fast.1way (lev, y, w = rep(1, length(y))) 
 
find.upcross (fun, fun.info, upcross.level = 0, guess = 1, tol = 
1e-05) 
 
gauss.cov (...) 


gcv.sreg (out, lambda.grid = NA, cost = 1, nstep.cv = 20, rmse = 
     NA, offset = 0, trmin = NA, trmax = NA, verbose = TRUE, tol = 1e-05, 
     find.min = TRUE, method = "GCV") 
 
golden.section.search (ax, bx, cx, f, niter = 25, f.extra = NA, 
     tol = 1e-05) 
 
grey.level (n) 
 
image.plot.info (...) 
image.plot.plt ( add=FALSE,  legend.shrink = 0.9, legend.width = 0.04, 
horizontal = FALSE, offset = 2 * legend.width, bigplot = NULL, smallplot = NULL)

image.smooth.setup (nrow = 64, ncol = 64, dx = 1, dy = 1, 
     kernel.function = function(x) {    exp(-abs(x))}, 
     theta = 1, Mwidth = nrow, Nwidth = ncol, ...) 
 
in.poly (xd, xp, convex.hull = FALSE, inflation=1e-7) 
 
krig.image.parameters (out) 
 
lonlat2xy (lnlt, miles = FALSE) 
 
make.surface.grid (grid.list, X, nx = 30, ny = 30, 
      info.list = FALSE, FUN = median) 
make.tmatrix (x, m = 2) 
make.tmatrix.F (x, m = 2) 
make.tmatrix.S (x, m = 2) 
 
matern (t , scale = 1, range = 1,alpha=1/range,
                smoothness = 0.5, nu= smoothness, phi=scale)  

minimax.crit (obj, des = TRUE, R) 
 
periodic.cov.1d (x1, x2, a, b) 
periodic.cov.cyl (x1, x2, a = 0, b = 365, theta = 1) 
periodic.plane.3d (x1, x2, a = 0, b = 365, theta = 1) 
 
plot.coef (x, cut.min = 8, graphics.reset = TRUE, common.range = FALSE) 
plot.krig.image (out, main = NA, digits = 4, which = rep(TRUE, 4), 
      graphics.reset = TRUE, ...) 
plot.qsreg (out, pch = "*", main = NA) 
plot.sim.krig.image (obj) 
plot.spatial.design (out) 
plot.vgram.matrix (obj, ...)
 
predict.interp.surface (obj, loc) 
predict.krig.image (obj, x, ...) 
predict.qsreg (out, x, derivative = 0, model = out$ind.cv.ps)
predict.sreg (out, x, derivative = 0, model = 1) 
 
print.krig.image (obj) 
print.qsreg (x, ...) 
print.spatial.design (x) 
print.sreg (x, ...) 
print.summary.Krig (x, ...) 
print.summary.krig.image (x, ...) 
print.summary.spatial.design (x, digits = 4) 
print.summary.sreg (x, ...) 
 
qr.q2ty (qr, y) 
 
qr.yq2 (qr, y) 
 
qsreg.fit (x, y, lam, maxit = 50, maxit.cv = 10, tol = 1e-04, 
    offset = 0, sc = sqrt(var(y)) * 1e-07, alpha = 0.5, wt = rep(1, 
    length(x)), cost = 1)
qsreg.psi ( r,alpha=.5,C=1) 
qsreg.rho ( r,alpha=.5,C=1)
qsreg.trace (x, y, lam, maxit = 50, maxit.cv = 10, tol = 1e-04, 
    offset = 0, sc = sqrt(var(y)) * 1e-07, alpha = 0.5, 
    wt = rep(1, length(x)), cost = 1) 

rad.cov (x1, x2, p = 1, with.log = TRUE, with.constant = TRUE, C =NULL) 
rad.cov.F (x1, x2, p = 1, with.log = TRUE, with.constant = TRUE, C = NULL) 
rad.cov.S (x1, x2, p = 1, with.log = TRUE, with.constant = TRUE) 
rad.covC (x1, x2, C, p = 1, with.log = TRUE, with.constant = TRUE) 
rad.image.cov (ind1, ind2, Y, cov.obj = NULL, setup = FALSE, grid,...) 
radbas.constant (m, d) 
 
rdist.F (x1, x2) 
rdist.S (x1, x2, lon.lat = FALSE) 
 
replace.args.function (fun, ...) 
 
sim.krig.image (out, nreps = 10) 
sim.rf.default (obj) 
 
sreg.df.to.lambda (df, x, wt, guess = 1, tol = 1e-05) 
sreg.fdf (h, info) 
sreg.fgcv (lam, obj) 
sreg.fgcv.model (lam, obj) 
sreg.fgcv.one (lam, obj) 
sreg.fit (lam, info) 
sreg.fs2hat (lam, obj) 
sreg.trace (h, info) 

stats.bplot (x, by, style = "tukey", outlier = TRUE) 
stats.sim.krig.image (obj) 
 
summary.gcv.Krig (out, lambda, cost = 1, verbose = FALSE, offset = 
    0, y = NULL) 
summary.gcv.sreg (out, lambda, cost = 1, nstep.cv = 20, offset = 
    0, verbose = TRUE) 
summary.krig.image (x, digits = 4, ...) 
summary.qsreg (x, ...) 
summary.spatial.design (x, digits = 4, ...) 
summary.sreg (x, digits = 4, ...) 

surface (x, grid.list = NA, extrap = FALSE, xlab = NULL, ylab = 
    ,zlab = NULL, main = NULL, levels = NULL, zlim = NULL, ...) 
surface.default (...) 
surface.krig.image (obj, grid.list = NA, extrap = TRUE, 
      graphics.reset = FALSE, xlab = NULL, ylab = NULL, main = NULL,
      zlab = NULL,zlim = NULL,levels = NULL, ptype = "I", ...) 
 
surface.surface (...)
 
test.data  

test.data2

unscale (x, x.center, x.scale) 
 
variogram.matrix (dat, R = 5, dx) 
 
world.dat 
}
\keyword{internal}
% docclass is function
% Converted by Sd2Rd version 1.21.
