## Criterion to minimize, proportional to the negative marginal log likelihood.
## Require (thetaTrain, XtX)
hML <- function(len) {
    l <- length(thetaTrain)
    k <- ncol(XtX) %/% l
    Jk <- J(k)
    message("parameter: ", len)
    K <- kerSEmat(thetaTrain, len = len)
    logdetK <- determinant(K, logarithm = TRUE)
    stopifnot(logdetK$sign == 1)
    Kinv <- solve(K)
    breveBox <- XtX * kronecker(Kinv, Jk) #7ms
    logdetBbox <- determinant(breveBox, logarithm = TRUE)
    stopifnot(logdetBbox$sign == 1)
    return(data.table(len = len, logdetK = logdetK$modulus, logdetBbox = logdetBbox$modulus))
}
## ## Prefers singular correlation matrix, i.e. infinite lengthscale.
## lenTest <- seq(0.1, 0.6, by = 0.1)
## system.time(DTll <- map_dfr(lenTest, hML)) #0.33s/6
## DTll[, plot(len, n * logdetK + logdetBbox, type = 'o')] #monotonic
## DTll[, plot(len, 60 * logdetK + logdetBbox, type = 'o')] #convex

