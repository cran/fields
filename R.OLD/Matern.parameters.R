Matern.cor.to.range <- function(d, nu, cor.target = 0.5, 
    guess = NULL, ...) {
    # define local function for root finding
    #
    ftemp <- function(theta, f.extra) {
        Matern(f.extra$d/theta, nu = f.extra$nu) - f.extra$cor.target
    }
    # inital guess is exponential
    if (is.null(guess)) {
        guess[1] <- guess[2] <- -d/log(cor.target)
    }
    #  extra info for function
    f.extra = list(d = d, nu = nu, cor.target = cor.target)
    # find  guesses that are above and below
    while (ftemp(guess[2], f.extra) < 0) {
        guess[2] <- guess[2] * 2
    }
    while (ftemp(guess[1], f.extra) > 0) {
        guess[1] <- guess[1]/2
    }
    temp <- bisection.search(guess[1], guess[2], f = ftemp, f.extra = f.extra, 
        ...)
    return(temp$x)
}
