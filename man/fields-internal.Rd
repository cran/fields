\name{fields internal }
\alias{fields internal}
\title{
  Fields internal and secondary functions 
}
\description{
Listed below are supporting fucntions for the major methods in fields. 
}
\usage{

COR (dat) 
 
D4.cov (m, n, x1 = m/2, x2 = n/2, weights = c(1, 0.1, 0.01, 
    0.001, 5e-04, 1e-04), cut.min = 8, max.m = m, max.n = n, 
    setup = F, reg.obj = NULL, grid = NULL) 
D4transform.image (x, inv = F, transpose = F, cut.min = 8) 
 
FIELDS.version  

Krig.coef (out, u = NULL, lambda = NULL, yM = NULL) 
Krig.df.to.lambda (df, D, guess = 1, tol = 1e-05) 
Krig.discretize (x, m = 64, n = 64, xg = NULL, yg = NULL, grid 
    = NULL,  expand = c(1, 1)) 
Krig.fdf (llam, info) 
Krig.fgcv (lam, obj) 
Krig.fgcv.model (lam, obj) 
Krig.fgcv.one (lam, obj) 
Krig.find.gcvmin (info, lambda.grid, gcv, gcv.fun, tol, verbose = F) 
Krig.fs2hat (lam, obj) 
Krig.ftrace (lam, D) 
Krig.newy (out, newy, verbose = F) 
Krig.replicates (out, verbose = F) 
Krig.updateY (out, Y, verbose = F) 
Krig.ynew (out, ynew, verbose = F) 
 
US (xlim = c(-124.7, -67.1), ylim = c(25.2, 49.4), add = F, ...) 
 
W.image.cov (ind1, ind2, Y, cov.obj = NULL, setup = F, m = 64, 
         n = 64, ...) 
 
WD4 (x) 
WD42d (x) 
WD42di (x) 
WD4i (x) 
 
WQS (x) 
WQS.T (x) 
WQS.basis (N, cut.n = 8) 
WQS2d (x, transpose = F) 
WQS2di (x, transpose = F) 
WQSdi (x, transpose = F) 
WQSi (x) 
WQSi.T (x) 
Wtransform.D (nx, ny, weights = c(1), cut.min = 8, details = F) 
 
bisection.search (x1, x2, f, tol = 1e-07, niter = 25, f.extra = 
        NA, upcross.level = 0) 
 
bplot.obj (data, pos = NA, width = NULL, labels = NULL, srt = 
            ,add = F, space = 0.25, sort.names = F, xlab = "", ylab = "", 
              label.cex = 1, xaxt = "n", outlier = T, horizontal = F, ...) 
 
cat.matrix (mat, digits = 8) 
 
cat.to.list (x, a) 
 
ceiling2 (m) 
 
conjugate.gradient (b, multAx, start = NULL, tol = 1e-05, kmax = 
             25, verbose = T,    ...) 
 
cover.design.S (R, nd, nruns = 1, nn = T, num.nn = 100, fixed = 
          ,scale.type = "range", R.center, R.scale, P = -20, Q = 20, 
          start = NULL, DIST = NULL, return.grid = T, return.transform = T) 
 
describe (x) 
describe.bplot (temp, style = "tukey", outlier = T) 
 
draw.bplot (temp, width, xpos, outlier = T, style = "tukey") 
draw.bplot.obj (obj, width, xpos, outlier = T, horizontal = F) 
 
dup (dat) 
 
dup.matrix (mat) 
 
dyn.load.fields (...) 
 
exp.cov.F (x1, x2, theta = rep(1, ncol(x1)), p = 1, C = NA) 
exp.cov.S (x1, x2, theta = 1) 
exp.earth.cov (x1, x2, theta = 1) 
 
fast.1way (lev, y, w = rep(1, length(y))) 
 
find.upcross (fun, fun.info, upcross.level = 0, guess = 1, tol = 
1e-05) 
 
gauss.cov (...) 
 
gcv.Krig (out, lambda.grid = NA, cost = 1, nstep.cv = 80, rmse = 
     NA,verbose = F, tol = 1e-05, offset = 0, y = NULL, lambda = NA) 
 
gcv.sreg (out, lambda.grid = NA, cost = 1, nstep.cv = 20, rmse = 
     NA, offset = 0, trmin = NA, trmax = NA, verbose = T, tol = 1e-05, 
     find.min = T, method = "GCV") 
 
golden.section.search (ax, bx, cx, f, niter = 25, f.extra = NA, 
     tol = 1e-05) 
 
grey.level (n) 
 
image.plot.info (...) 
image.plot.plt ( add=F,  legend.shrink = 0.9, legend.width = 0.04, 
horizontal = F, offset = 2 * legend.width, bigplot = NULL, smallplot = NULL)

image.smooth.setup (nrow = 64, ncol = 64, dx = 1, dy = 1, 
     kernel.function = function(x) {    exp(-abs(x))}, 
     theta = 1, Mwidth = nrow, Nwidth = ncol, ...) 
 
in.poly (xd, xp, convex.hull = F) 
 
krig.image.parameters (out) 
 
lonlat2xy (lnlt, miles = F) 
 
make.Amatrix.Krig (out, x0 = out$x, lambda) 
make.surface.grid (grid.list, X, nx = 30, ny = 30, 
      info.list = F, FUN = median) 
make.tmatrix (x, m = 2) 
make.tmatrix.F (x, m = 2) 
make.tmatrix.S (x, m = 2) 
 
matern (x = seq(0, 4 * range, len = 100), scale = 1, range = 1, 
    smoothness = 0.5) 
 
minimax.crit (obj, des = T, R) 
 
nkden (data, bandwidth, n.points, grid) 
 
periodic.cov.1d (x1, x2, a, b) 
periodic.cov.cyl (x1, x2, a = 0, b = 365, theta = 1) 
periodic.plane.3d (x1, x2, a = 0, b = 365, theta = 1) 
 
plot.coef (x, cut.min = 8, graphics.reset = T, common.range = F) 
plot.krig.image (out, main = NA, digits = 4, which = rep(T, 4), 
      graphics.reset = T, ...) 
plot.qsreg (out, pch = "*", main = NA) 
plot.sim.krig.image (obj) 
plot.spatial.design (out) 
 
predict.interp.surface (obj, loc) 
predict.krig.image (obj, x, ...) 
predict.qsreg (out, x, derivative = 0, model = 1) 
predict.sreg (out, x, derivative = 0, model = 1) 
predict.tps (out, x, y, lambda, df, omega, u = NA, 
      derivative = 0, model = NA) 
 
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
 
qsreg.trace (x, y, lam, maxit = 50, maxit.cv = 10, tol = 1e-04, 
    offset = 0, sc = sqrt(var(y)) * 1e-07, alpha = 0.5, 
    wt = rep(1, length(x)), cost = 1) 
 
rad.cov (x1, x2, p = 1, with.log = T, with.constant = T, C =NULL) 
rad.cov.F (x1, x2, p = 1, with.log = T, with.constant = T, C = NULL) 
rad.cov.S (x1, x2, p = 1, with.log = T, with.constant = T) 
rad.covC (x1, x2, C, p = 1, with.log = T, with.constant = T) 
rad.image.cov (ind1, ind2, Y, cov.obj = NULL, setup = F, grid,...) 
radbas.constant (m, d) 
 
rdist.F (x1, x2) 
rdist.S (x1, x2, lon.lat = F) 
 
replace.args.function (fun, ...) 
 
setup.W.image.cov (m, n, weights = 1, D = NULL, cut.min = 8, 
     max.m = m,max.n = n, x = NULL, grid = NULL, normalize = T, ...) 
 
setup.Wtransform.cov (m, n, weights = 1, cut.min = 8, max.m = m, 
     max.n = n,x = NULL, grid = NULL, normalize = T, ...) 
 
sim.krig.image (out, nreps = 10) 
sim.rf.W (cov.obj) 
sim.rf.default (obj) 
 
sreg.df.to.lambda (df, x, wt, guess = 1, tol = 1e-05) 
sreg.fdf (h, info) 
sreg.fgcv (lam, obj) 
sreg.fgcv.model (lam, obj) 
sreg.fgcv.one (lam, obj) 
sreg.fit (lam, info) 
sreg.fs2hat (lam, obj) 
sreg.trace (h, info) 

stats.bplot (x, by, style = "tukey", outlier = T) 
stats.sim.krig.image (obj) 
 
summary.Krig (x, digits = 4, ...) 
summary.gcv.Krig (out, lambda, cost = 1, verbose = F, offset = 
    0, y = NULL) 
summary.gcv.sreg (out, lambda, cost = 1, nstep.cv = 20, offset = 
    0, verbose = T) 
summary.krig.image (x, digits = 4, ...) 
summary.qsreg (x, ...) 
summary.spatial.design (x, digits = 4, ...) 
summary.sreg (x, digits = 4, ...) 

surface (x, grid.list = NA, extrap = F, xlab = NULL, ylab = 
    ,zlab = NULL, main = NULL, levels = NULL, zlim = NULL, ...) 
surface.default (...) 
surface.krig.image (obj, grid.list = NA, extrap = T, 
      graphics.reset = F, xlab = NULL, ylab = NULL, main = NULL,
      zlab = NULL,zlim = NULL,levels = NULL, ptype = "I", ...) 
 
surface.surface (obj, lab = NA, type = "b", zlab, xlab, ylab, 
   graphics.reset = F,  ...) 
 
test.data  

test.data2

unique.matrix (x) 
 
unscale (x, x.center, x.scale) 
 
variogram.matrix (dat, R = 5) 
 
world.dat 
}
\keyword{misc}
% docclass is function
% Converted by Sd2Rd version 1.21.
