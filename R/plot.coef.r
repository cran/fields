"plot.coef" <-
function (x, cut.min = 8, graphics.reset = T, common.range = F) 
{
    old.par <- par()
    par(mar = c(3, 0, 3, 0))
    n <- dim(x)[1]
    m <- dim(x)[2]
    NN <- n
    MM <- m
    level <- 1
    while (min(c(NN, MM)) > cut.min) {
        level <- level + 1
        NN <- NN/2
        MM <- MM/2
    }
    print(level)
    print(MM)
    print(NN)
    n2 <- NN
    m2 <- MM
    n3 <- NN * 2
    m3 <- MM * 2
    n1 <- 1
    m1 <- 1
    zr <- range(x)
    set.panel(level, 3)
    image.plot(x[m1:m2, n1:n2], zlim = zr, xaxt = "n", yaxt = "n", 
        graphics.reset = F)
    mtext(3, line = 1, text = "smooth", cex = 1.1)
    par(mfg = c(2, 1, level, 3))
    level <- 1
    while (n3 <= n & m3 <= m) {
        cat(c(m1, m2, m3), fill = T)
        cat(c(n1, n2, n3), fill = T)
        zr <- range(c(x[m1:m2, (n2 + 1):n3], x[(m2 + 1):m3, n1:n2], 
            x[(m2 + 1):m3, (n2 + 1):n3]))
        image(x[m1:m2, (n2 + 1):n3], zlim = zr, xaxt = "n", yaxt = "n")
        mtext(2, line = 1, text = paste("level", level, "detail"), 
            cex = 1.1)
        if (level == 1) {
            mtext(3, line = 1, text = "Horizontal", cex = 1.1)
        }
        image(x[(m2 + 1):m3, n1:n2], zlim = zr, xaxt = "n", yaxt = "n")
        if (level == 1) {
            mtext(3, line = 1, text = "Vertical", cex = 1.1)
        }
        image.plot(x[(m2 + 1):m3, (n2 + 1):n3], zlim = zr, xaxt = "n", 
            yaxt = "n")
        if (level == 1) {
            mtext(3, line = 1, text = "Diagonal", cex = 1.1)
        }
        level <- level + 1
        n2 <- n3
        n3 <- n2 * 2
        m2 <- m3
        m3 <- m2 * 2
    }
    if (graphics.reset) {
        par(old.par)
    }
    invisible()
}
