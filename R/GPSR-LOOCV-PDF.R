#' Log determinant of predictive "covariance matrix"
#' \eqn{log(det(\tilde{\Sigma}))}
LogDetCov <- function(n, sigma2, eps2) {
    sum(log(sigma2 + eps2)) + (n - length(sigma2)) * log(eps2)
}

#' Log determinant of "inverse match matrix"
#' \eqn{log(det(X^T \tilde{\Sigma}^{-1} X))}
LogDetXSinvX <- function(X, V, sigma2, eps2) {
    n <- ncol(V)
    VtX <- Matrix::crossprod(V, X)
    XSinvX <- Matrix::crossprod(VtX, Diagonal(x = 1 / (sigma2 + eps2)) %*% VtX)
    R <- X - V %*% VtX
    XSinvX <- XSinvX + Matrix::crossprod(R) / eps2
    determinant(XSinvX, logarithm = TRUE)
}

## Compute LOO log determinants from scratch: "covariance", "inverse match"
FunPred <- function(i, len) {
    l <- length(thetaTrainAll) - 1
    thetaTrain <- thetaTrainAll[-i]
    Xdi <- Matrix::Matrix(AXtrain[,,-i], n, k * l)
    force(len)
    ret <- GPSubspacePrep(thetaTrain, Xdi, len) #1.22s
    list2env(ret, environment(GPSubspacePred))
    ## Require (K, Kinv, XtX, svdX) in .GlobalEnv
    ret <- GPSubspacePred(thetaTrainAll[i], thetaTrain, len, NULL) #1.2s (k=20,l=10)
    ## str(ret)
    ## plot(ret$sigma2, log = 'y')
    ## plot(ret$sigma2[seq(30)], log = 'y')
    ## abline(h = ret$eps2)
    ## abline(v = k)
    ## plot(log(ret$sigma2 + ret$eps2))
    denom <- LogDetCov(n, ret$sigma2, ret$eps2) #-5.3e5
    ## LogDetCov(n, ret$sigma2, ret$eps2) - n * log(ret$eps2) #340
    enum <- LogDetXSinvX(AXtrain[,,i], ret$V, ret$sigma2, ret$eps2)
    ## 1 - -logdetcrossBox[i] / enum$modulus
    return(data.table(theta = thetaTrainAll[i], len = len, enum = enum$modulus, denom = denom))
}
## ## NOTE: Log LOOCV prob are numerically the same as those computed via `GPSubspacePred`,
## ## LOOCV prob prefers large lengthscale:
## ## (1) the "covariance matrix" shrinks in volume, which boosts prob density;
## ## (2) the "mismatch matrix" grows, which reduces prob density, but is less influential than (1).
## ## Log LOOCV prob can still be useful if n is not way larger than k.
## lenTest <- 0.25
## DTpred <- map_dfr(seq(thetaTrainAll), FunPred, len = lenTest)
## DT[len == lenTest, h] / DTpred[, sum(denom * k + enum * n) / 2] - 1
## ## Matched log LOOCV prob: len = 0.1, 0.15, 0.2, 0.25, 0.3, 0.35
## DTscan <- map_dfr(seq(0.1, 0.5, by = 0.1), ~FunPred(i = 6, len = .))
## DTscan[, denom * k + enum * n]

#' Log determinant of leave-one-out predictive "covariance matrix"
#' \eqn{log(det(\tilde{\Sigma}_{-i}))}
#' Result does not match raw computation!
#' @note require (X, XtX; Kinv)
LogDetLooCov <- function(i) {
    k <- ncol(XtX) %/% ncol(Kinv)
    Jk <- J(k)
    kbii <- Kinv[i, i]
    kbio <- Kinv[-i, i]
    Kbi <- Kinv[-i, -i]
    ## Deltai <- Kbi / kbio / rep(kbio, each = l-1) * kbii - 1
    Dkbioi <- Diagonal(x = 1 / kbio)
    Deltai <- Dkbioi %*% Kbi %*% Dkbioi * kbii - 1
    ## This can fail.
    ## Deltai <- as(Deltai, "dsyMatrix")
    ## This takes a trianglar submatrix and call it symmetric.
    Deltai <- forceSymmetric(Deltai)
    ## Row/colum indices of the i-th training point.
    idxi <- (i-1)*k + seq(k)
    deltaBoxi <- XtX[-idxi, -idxi] * kronecker(Deltai, Jk)
    deltaBoxi <- as(deltaBoxi, "dpoMatrix")
    ## Note: Sigma can be explicitly computed as follows, but it takes a lot of time!
    ## Sigma <- (Diagonal(n) + X[,-idxi] %*% tcrossprod(deltaBoxiInv, X[,-idxi])) / kbii
    U <- chol(deltaBoxi)
    Uinv <- solve(U)
    XiUinv <- X[,-idxi] %*% Uinv #0.154s
    ret <- svd(XiUinv, nu = 0, nv = 0) #0.326s
    sum(log(ret$d^2 + 1)) - n * log(kbii)
}

#' Log determinant of the principal submatrices of breveBoxinv
#'
#' \eqn{log(det(crossBox_i)) = -log(det(X_i^T \tilde{\Sigma}_{-i}^{-1} X_i))}
#' @param breveBoxinv inverse of breveBox
#' @return a vector of length l, one for each \eqn{log(det(crossBox_i))}.
#' @note Require (k)
LogDetCrossBox <- function(breveBoxinv) {
    l <- ncol(breveBoxinv) %/% k
    ret <- double(l)
    for(i in seq(l)) {
        idxi <- (i-1)*k + seq(k)
        detcrossBoxi <- determinant(breveBoxinv[idxi, idxi], logarithm=TRUE) #5ms
        ## stopifnot(detcrossBoxi$sign == 1)
        if(detcrossBoxi$sign == -1) warning("detcrossBoxi is negative: i = ", i)
        ret[i] <- detcrossBoxi$modulus
    }
    ret
}

#' Leave-one-out predictive probability density
#' Faster than using `LogDetLooCov()` and `LogDetCrossBox()`.
#' @param i index of left out parameter point
#' @param Kinv inverse of correlation matrix
#' @param logdetbreveBox log det of breveBox
#' @param logdetcrossBoxi log det of crossBox_i, computed from `LogDetCrossBox()`.
#' @note Require global (XtX)
LogLooPredProb <- function(i, Kinv, logdetbreveBox, logdetcrossBoxi) {
    l <- ncol(Kinv)
    k <- ncol(XtX) %/% l
    Jk <- J(k)
    kbii <- Kinv[i, i]
    kbio <- Kinv[-i, i]
    Kbi <- Kinv[-i, -i]
    vi <- -kbio / kbii
    Dkbioi <- Diagonal(x = 1 / kbio)
    Deltai <- Dkbioi %*% Kbi %*% Dkbioi * kbii - 1
    Deltai <- forceSymmetric(Deltai)

    ## Row/colum indices of the i-th training point.
    idxi <- (i-1)*k + seq(k)
    deltaBoxi <- XtX[-idxi, -idxi] * kronecker(Deltai, Jk) #8ms
    message("Converting deltaBoxi: i = ", i)
    deltaBoxi <- as(deltaBoxi, "dpoMatrix") #4ms
    detdeltaBoxi <- determinant(deltaBoxi, logarithm=TRUE)
    ## stopifnot(detdeltaBoxi$sign == 1)
    if(detdeltaBoxi$sign == -1) warning("detdeltaBoxi is negative: i = ", i)
    logdetdeltaBoxi <- detdeltaBoxi$modulus
    ## logdetcrossBoxi <- LogDetCrossBox(i)
    ## Putting everything together.
    logdetLooCovi <- (- (n + k * (l-1)) * log(kbii)
        - 2 * k * log(prod(abs(vi)))
        + (logdetbreveBox - logdetdeltaBoxi)
        + logdetcrossBoxi)
    ## return(logdetLooCovi)
    return(-0.5 * (k * logdetLooCovi - n * logdetcrossBoxi))
    ## return(-0.5 * (
    ##     - k * (n + k * (l-1)) * log(kbii)
    ##     - 2 * k^2 * log(prod(abs(vi)))
    ##     + k * (logdetbreveBox - logdetdeltaBoxi)
    ##     + (k - n) * logdetcrossBoxi
    ## ))
}
## system.time(logdetTildeSigma <- map_dbl(seq(l), LogDetLooCov)) #3.6s
## system.time(logdetcrossBox <- LogDetCrossBox(breveBoxinv)) #0.024s
## logpiSlow <- -0.5 * (k * logdetTildeSigma - n * logdetcrossBox)
## system.time(logpi <- map2_dbl(seq(l), logdetcrossBox,
##                               ~LogLooPredProb(.x, Kinv, logdetbreveBox, .y))) #2.1s
## all.equal(logpi, logpiSlow)

## Compute LOO log determinants: "covariance", "inverse match"
FunLOO <- function(thetaTrain, XtX, len) {
    l <- length(thetaTrain)
    K <- kerSEmat(thetaTrain, len = len)
    Kinv <- solve(K)
    breveBox <- XtX * kronecker(Kinv, Jk) #7ms
    breveBoxinv <- solve(breveBox) #2ms
    logdetcrossBox <- LogDetCrossBox(breveBoxinv) #0.024s
    detbreveBox <- determinant(breveBox, logarithm = TRUE) #2ms
    if(detbreveBox$sign == -1) warning("detbreveBox is negative")
    logdetbreveBox <- detbreveBox$modulus
    logdetLooCov <- map2_dbl(seq(l), logdetcrossBox, ~LogLooPredProb(.x, Kinv, logdetbreveBox, .y))
    data.table(theta = thetaTrain,
               logdetLooXSinvX = -logdetcrossBox,
               logdetLooCov = logdetLooCov)
}
## DTloo <- FunLOO(thetaTrainAll, Matrix::crossprod(X), len)
## DTloo[, logdetLooXSinvX * n + logdetLooCov * k] / DTpred[, denom * k + enum * n] - 1
## DTloo[, plot(theta, logdetLooXSinvX * n + logdetLooCov * k)]


#' Criterion to minimize, the negative log LOO predictive probability.
#' @note require global (thetaTrain, XtX)
h <- function(len) {
    l <- length(thetaTrain)
    k <- ncol(XtX) %/% l
    Jk <- J(k)
    message("parameter: ", len)
    K <- kerSEmat(thetaTrain, len = len)
    Kinv <- solve(K)
    breveBox <- XtX * kronecker(Kinv, Jk) #7ms
    ## message("Converting breveBox")
    ## breveBox <- as(breveBox, "dpoMatrix") #1ms
    breveBoxinv <- solve(breveBox) #2ms
    detbreveBox <- determinant(breveBox, logarithm = TRUE) #2ms
    ## stopifnot(detbreveBox$sign == 1)
    if(detbreveBox$sign == -1) warning("detbreveBox is negative")
    logdetbreveBox <- detbreveBox$modulus
    logdetcrossBox <- LogDetCrossBox(breveBoxinv) #0.024s
    logpi <- map2_dbl(seq(l), logdetcrossBox, ~LogLooPredProb(.x, Kinv, logdetbreveBox, .y)) #2.1s
    -sum(logpi)
}
## system.time(ret <- h(0.5)) #2.1s
## system.time(ret <- optim(0.2, h, method = "L-BFGS-B", lower = 0.1, upper = 3)) #10s
## ## system.time(ret <- optim(1, h, method = "L-BFGS-B", lower = 1, upper = 3)) #10s
## system.time(ret <- optimize(h, c(0.1, 0.6)))
## system.time(ret <- optimize(h, c(0.1, 0.35)))
## system.time(ret <- optimize(h, c(0.34, 0.39)))
## system.time(ret <- optimize(h, c(0.387, 0.4)))
## system.time(ret <- optimize(h, c(0.39, 0.43)))
## DT <- data.table(len = seq(0.1, 0.35, by = 0.05))
## DT <- data.table(len = seq(0.4, 0.5, by = 0.05))
## DT$h <- NA_real_
## for(i in seq(nrow(DT))) {
##     DT[i, h := h(len)]
##     ## DT[i, h2 := h2(len)]
## }
