"cover.design" <-
function (R, nd, nruns = 1, nn = TRUE, num.nn = 100, fixed = NULL, 
    scale.type = "unscaled", R.center, R.scale, P = -20, Q = 20, 
    start = NULL, DIST = NULL, return.grid = TRUE, return.transform = 
TRUE,max.loop=20) 
{
    if (!is.null(start) && is.matrix(start)) {
        if (any(duplicated.array(start))) 
            stop("Error: start must not have duplicate rows")
        start <- rowmatch(start, R)
        if (any(is.na(start))) 
            stop("Error: Starting design must be a subset of R")
    }
    R.orig <- R
    R <- as.matrix(R)
    if (any(duplicated.array(R))) 
        stop("Error: R must not have duplicate rows")
    if (num.nn >= nrow(R)) 
        nn <- FALSE
    if (is.null(DIST)) 
        DIST <- function(x, y) {
            rdist(x, y)
        }
    id <- 1:nrow(R)
    if (!is.null(start)) 
        nd <- length(start)

    if (is.null(fixed)) 
        n <- nd
    else {
        n <- nd + length(fixed)
    }

    R <- transformx(R, scale.type, R.center, R.scale)
    transform <- attributes(R)
    saved.crit <- rep(NA, nruns)
    saved.designs <- matrix(NA, nrow = nruns, ncol = n)
    saved.hist <- list(1:nruns)
    for (RUNS in 1:nruns) {
        if (is.null(start)) {
            if (!is.null(fixed)) {
                Dset <- sample((1:nrow(R))[-fixed], nd)
                Dset <- c(Dset, fixed)
            }
            else Dset <- sample(1:nrow(R), nd)
        }
        else {
            if (length(start) > nd) 
                stop("Error: the start matrix must have nd rows")
            Dset <- start
            if (!is.null(fixed)) 
                Dset <- c(Dset, fixed)
        }
        design.original <- R.orig[Dset, ]
        Dset.orginal <- Dset
        Cset <- id[-Dset]
        dist.mat <- DIST(R[Cset, ], R[Dset, ])
        rs <- dist.mat^P %*% rep(1, n)
        crit.i <- crit.original <- sum(rs^(Q/P))^(1/Q)
        CRIT <- rep(NA, length(Cset))
        CRIT.temp <- rep(NA, length(Cset))
        hist <- matrix(c(0, 0, crit.i), ncol = 3, nrow = 1)
        loop.counter<-1
        repeat {
        for (i in 1:nd) {
                Dset.i <- matrix(R[Dset[i], ], nrow = 1)
                partial.newrow <- sum(DIST(Dset.i, R[Dset[-i], 
                  ])^P)
                rs.without.i <- rs - c(DIST(Dset.i, R[-Dset, 
                  ])^P)
                if (nn) 
                  vec <- (1:length(Cset))[order(dist.mat[, i])[1:num.nn]]
                else vec <- 1:length(Cset)
                for (j in vec) {
                  Cset.j <- matrix(R[Cset[j], ], nrow = 1)
                  newcol <- c(DIST(Cset.j, R[c(-Dset, -Cset[j]), 
                    ])^P)
                  CRIT[j] <- (sum((rs.without.i[-j] + newcol)^(Q/P)) + 
                    (DIST(Cset.j, Dset.i)^P + partial.newrow)^(Q/P))^(1/Q)
                }
                best <- min(CRIT[!is.na(CRIT)])
                best.spot <- Cset[CRIT == best][!is.na(Cset[CRIT == 
                  best])][1]
                crit.old <- crit.i
                if (best < crit.i) {
                  crit.i <- best
                  hist <- rbind(hist, c(Dset[i], best.spot, crit.i))
                  Dset[i] <- best.spot
                  Cset <- id[-Dset]
                  dist.mat <- DIST(R[Cset, ], R[Dset, ])
                  rs <- (dist.mat^P) %*% rep(1, n)
                }
            }
            if ((crit.i == crit.old) |(loop.counter>=max.loop))
                break
        loop.counter<- loop.counter +1
       
        }
        saved.crit[RUNS] <- crit.i
        saved.designs[RUNS, ] <- Dset
        saved.hist[[RUNS]] <- hist
    }
    ret <- (1:nruns)[saved.crit == min(saved.crit)]
    if (length(ret) > 1) {
        print("Greater than 1 optimal design; keeping first one......")
        ret <- ret[1]
    }
    crit.i <- saved.crit[ret]
    hist <- saved.hist[[ret]]
    nhist <- nrow(hist)
    nloop <- nruns
    hist <- cbind(c(0:(nrow(hist) - 1)), hist)
    dimnames(hist) <- list(NULL, c("step", "swap.out", "swap.in", 
        "new.crit"))
    out.des <- R[saved.designs[ret, ], ]
    out.des <- unscale(out.des, transform$x.center, transform$x.scale)
    out <- list(design = out.des, call = match.call(), best.id = c(saved.designs[ret, 
        ]), fixed = fixed, opt.crit = crit.i, start.design = design.original, 
        start.crit = crit.original, history = hist, other.designs = saved.designs, 
        other.crit = saved.crit, DIST = DIST, nn = nn, num.nn = num.nn, 
        P = P, Q = Q, nhist = nhist - 1, nloop = (nloop - 1)/n)
    if (return.grid) 
        out$grid <- R.orig
    if (return.transform) 
        out$transform <- transform
    class(out) <- "spatial.design"
    out
}
