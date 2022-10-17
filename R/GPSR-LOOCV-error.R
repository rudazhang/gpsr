## Compute LOO Riemannian distances from scratch
## Require: (thetaTrain, AX; k)
FunDist <- function(i, len) {
    l <- length(thetaTrain) - 1
    thetaLOO <- thetaTrain[-i]
    Xdi <- Matrix::Matrix(AX[,,-i], n, k * l)
    force(len)
    ret <- GPSubspacePrepSVD(Xdi) #1.22s
    ## list2env(ret, environment(GPSubspacePredEVD))
    list2env(ret, .GlobalEnv)
    ## Require (K, Kinv, XtX, svdX) in .GlobalEnv
    K <- kerSEmat(thetaLOO, len = len)
    Kinv <- solve(K)
    list2env(mget(c("K", "Kinv")), .GlobalEnv)
    ret <- GPSubspacePredEVD(thetaTrain[i], thetaLOO, len, t=NULL) #1.2s (k=20,l=10)
    ret <- svd(Matrix::crossprod(ret$V[,seq(k)], AX[,,i]))
    sigma <- ret$d
    sigma[sigma > 1] <- 1
    angle <- acos(sigma)
    return(data.table(theta = thetaTrain[i], len = len, dist = vecnorm(angle)))
}
## lenTest <- seq(0.1, 0.5, by = 0.1)
## DTscan <- CJ(len = lenTest, i = seq(thetaTrain))
## system.time(DTdist <- map2_dfr(DTscan$i, DTscan$len, FunDist))
## ## 107s/5 (k=20,l=10); 209s/5(k=40,l=11)
## ## lenTest2 <- seq(0.22, 0.24, by = 0.02)
## ## lenTest2 <- c(0.35)
## ## DTscan2 <- CJ(len = lenTest2, i = seq(thetaTrain))
## ## system.time(DTdist2 <- map2_dfr(DTscan2$i, DTscan2$len, FunDist)) #83s/2
## DTstat <- DTdist[, .(d2 = vecnorm(dist), d1 = sum(dist)), keyby = .(len)]
## ## DTstat <- rbind(DTstat, DTdist2[, .(d2 = vecnorm(dist), d1 = sum(dist)), keyby = .(len)])
## setkey(DTstat, len)
## ## k=20,l=10: Both 1-norm and 2-norm approximately selects len=0.19.
## ## k=40,l=11: Both 1-norm and 2-norm approximately selects len=0.25.
## DTstat[, plot(len, d2 / max(d2), ylim = c(0.9, 1), type = 'o', col = "red")]
## DTstat[, lines(len, d1 / max(d1), type = 'o')]
## ## fwrite(DTstat, "LOOCV-dist-k40.csv")

## Compute LOO Riemannian distances: SVD version
## Require: (XtX)
LOODistSvd <- function(i, Kinv, VbtX) {
    k <- ncol(XtX) %/% ncol(Kinv)
    Jk <- J(k)
    ## Construct \Pi_{-1}
    kbii <- Kinv[i, i]
    kbio <- Kinv[-i, i]
    Kbi <- Kinv[-i, -i]
    Dkbioi <- Diagonal(x = 1 / kbio)
    Deltai <- Dkbioi %*% Kbi %*% Dkbioi * kbii - 1
    Deltai <- forceSymmetric(Deltai)
    ## Row/colum indices of the i-th training point.
    idxi <- (i-1)*k + seq(k)
    deltaBoxi <- XtX[-idxi, -idxi] * kronecker(Deltai, Jk) #0.02s
    deltaBoxi <- as(deltaBoxi, "dpoMatrix") ## 0.01s, computes Cholesky factor
    ## Triangular inverse
    U <- chol(deltaBoxi) #0ms
    Uinv <- solve(U) #4ms
    VXiUinv <- VbtX[,-idxi] %*% Uinv #3ms
    Vcheck <- irlba::irlba(VXiUinv, k) #0.11s (full svd returning k u's: 0.132s)
    sigma <- svd(Matrix::crossprod(VbtX[,idxi], Vcheck$u))$d #5ms
    sigma[sigma > 1] <- 1
    angle <- acos(sigma)
    vecnorm(angle)
}

#' Compute LOO prediction: EVD version
#' @return LOO prediction in the global basis, \equation{\circ{V}_{-i}}.
#' @note Require: (XtX; k, Jk)
LOOPredEvd <- function(i, Kinv, VbtX) {
    ## Row/colum indices of the i-th training point.
    idxi <- (i-1)*k + seq(k)
    ## Construct \Delta_{-i}, positive semi-definite
    kbii <- Kinv[i, i]
    kbio <- Kinv[-i, i]
    Kbi <- Kinv[-i, -i]
    Dkbioi <- Diagonal(x = 1 / kbio)
    Deltai <- Dkbioi %*% Kbi %*% Dkbioi * kbii - 1
    Deltai <- forceSymmetric(Deltai)
    ## Construct \Pi_{-1}, positive semi-definite
    PIi <- XtX[-idxi, -idxi] * kronecker(Deltai, Jk) #0.02s
    ## Add a nugget for numerical stability
    maxPIi <- max(abs(diag(PIi)))
    PIi <- PIi + Matrix::Diagonal(k * (l-1), epsMachine * maxPIi) #4ms
    PIi <- as(PIi, "dpoMatrix") ## 0.01s, computes Cholesky factor
    ## Construct P_i = A_{-i} (\Pi_{-i})^{-1} A_{-i}^T
    Ami <- VbtX[,-idxi]
    Sol <- solve(PIi, t(Ami)) #0.01s
    Pi <- Ami %*% Sol #3ms

    ## Added for the special case of G_{1,2}.
    if (nrow(Pi) == 2) {
        evd <- eigen(Pi, symmetric = TRUE)
        return(evd$vectors[, seq(k)])
    } else {
        evd <- RSpectra::eigs_sym(Pi, k) #8ms (full EVD, 88ms)
        return (evd$vectors)
    }
}

#' Compute LOO Riemannian distances: EVD version
#' @note Require: (XtX; k, Jk)
LOODistEvd <- function(i, Kinv, VbtX) {
    Vcirc <- LOOPredEvd(i, Kinv, VbtX)
    ## Row/colum indices of the i-th training point.
    idxi <- (i-1)*k + seq(k)

    ## sigma <- svd(Matrix::crossprod(VbtX[,idxi], Vcirc))$d #5ms
    ## sigma[sigma > 1] <- 1
    ## angle <- acos(sigma)
    ## vecnorm(angle)
    RiemannianDist(VbtX[,idxi], Vcirc)
}

#' Gradient of squared LOOCV Riemannian distance
#' @param lrGKinv a list of order-l matrices that consists of \eqn{\partial k_{p,q} / k_{p,q}},
#' if using isotropic lengthscales; or a list of them, if using separable lengthscales.
#' @param trunk Extended truncated EVD for approximate pseudo-inverse, defaults to 2*k, max r.
#' @return a vector for separable lengthscales, or a scalar for isotropic lengthscales.
#' @note Require: (XtX; k, l, Jk)
#' @details If `trunk` is close to r, `RSpectra::eigs_sym()` throws a warning that `eigen()` is used.
gradLOODist2Evd <- function(i, Kinv, lrGKinv, VbtX, trunk) {
    ## Duplicate part with `LOOPredEvd()`
    ## Row/colum indices of the i-th training point.
    idxi <- (i-1)*k + seq(k)
    ## Construct \Delta_{-i}, positive semi-definite
    kbii <- Kinv[i, i]
    kbio <- Kinv[-i, i]
    Kbi <- Kinv[-i, -i]
    Dkbioi <- Diagonal(x = 1 / kbio)
    Deltai <- Dkbioi %*% Kbi %*% Dkbioi * kbii - 1
    Deltai <- forceSymmetric(Deltai)
    ## Construct \Pi_{-1}, positive semi-definite
    XtXi <- XtX[-idxi, -idxi] #15ms (Suprise! This actually takes quite some time.)
    PIi <- XtXi * kronecker(Deltai, Jk) #13ms (23s in one line)
    PIi <- as(PIi, "dpoMatrix") ## 0.01s, computes Cholesky factor
    ## Construct P_i = A_{-i} (\Pi_{-i})^{-1} A_{-i}^T
    Ami <- VbtX[,-idxi]
    Sol <- solve(PIi, t(Ami)) #0.01s
    Pi <- Ami %*% Sol #3ms

    ## Extended truncated eigen-decomposition for approximate pseudo-inverse
    evd <- RSpectra::eigs_sym(Pi, trunk) #30ms (k EVD, 10ms; full EVD, 88ms)
    eigs <- evd$values
    Vcheck <- evd$vectors
    Vk <- Vcheck[, seq(k)]
    evp <- tcrossprod(rep(1, trunk), eigs[seq(k)])
    evShiftRev <- 1 / (evp - eigs)
    diag(evShiftRev) <- 0
    ## SVD for principal angles
    Ai <- VbtX[,idxi]
    svd <- svd(Matrix::crossprod(Ai, Vk)) #5ms
    sigma <- svd$d
    sigma[sigma > 1] <- 1
    angle <- acos(sigma)
    ## Partial gradient of squared LOOCV Riemannian distance
    ## Direct form may cause large error, ues series form for sigma = 1 - x, x < 1e-5.
    d <- 1 - sigma
    gd2s <- ifelse(d > 1e-5, angle / sqrt(1 - sigma^2),
                   1 + d/3 + 2/15 * d^2 + 2/35 * d^3 + 8/315 * d^4 + 8/693 * d^5)

    ## Partial derivative of squared LOOCV Riemannian distance w.r.t. one lengthscale
    ## Require: Deltai, XtXi, Sol; evd (eigs, Vcheck, Vk), evp (evShiftRev), svd (u, v), gd2s;
    ## i, l, k, Jk.
    pdd2 <- function(rGKinv) {
        ## Gradient of \Delta_{-i}, not positive semi-definite
        rGKinvVi <- tcrossprod(rGKinv[-i, i], rep(1, l-1))
        rGKinvMi <- rGKinv[-i, -i] + rGKinv[i, i] - rGKinvVi - t(rGKinvVi)
        rGKinvMi <- Matrix::forceSymmetric(rGKinvMi)
        GDeltai <- (Deltai + 1) * rGKinvMi
        ## Gradient of \Pi_{-i}, not positive semi-definite
        GPIi <- XtXi * kronecker(GDeltai, Jk) #13ms
        ## Gradient of top-k eigenvectors: 7ms (for loop, 123ms)
        GPiVk <- Matrix::crossprod(-Sol, GPIi %*% (Sol %*% Vk)) #4ms; v
        VtGPiVk <- Matrix::crossprod(Vcheck, GPiVk) #u
        wVtGPiVk <- (evShiftRev - 1 / evp) * VtGPiVk #w
        GVk <- Vcheck %*% wVtGPiVk + GPiVk %*% Diagonal(x = 1 / eigs[seq(k)])
        ## Gradient of singular values: 6ms (for loop, 11ms)
        GSV <- function(.) as.vector(Matrix::crossprod(Ai %*% svd$u[, .], GVk %*% svd$v[, .]))
        Gsigma <- vapply(seq(k), GSV, double(1))
        ## Gradient of squared LOOCV Riemannian distance
        -2 * sum(gd2s * Gsigma)
    }
    gd2 <- vapply(lrGKinv, pdd2, double(1))
    gd2
}

#' Criterion to minimize, the LOOCV prediction error: sum of squared Riemannian distances.
#' @note require (thetaTrain, XtX, VbtX; getKinv)
#' @export
hSSDist <- function(len) {
    message(paste(c("[objective]", "parameter: ", format(len, digits = 8)), collapse = "\t"))
    Kinv <- getKinv(thetaTrain, len)
    l <- ifelse(is.matrix(thetaTrain), nrow(thetaTrain), length(thetaTrain))
    dist <- vapply(seq_len(l), function(.) LOODistEvd(., Kinv, VbtX), double(1))
    sum(dist^2)
}
## system.time(ret <- hSSDist(0.5)) #0.826s

#' Gradient of sum of squared LOOCV Riemannian distances
#' @param len Correlation lengthscales.
#' @param trunk Extended truncated EVD for approximate pseudo-inverse, defaults to 2*k, max r.
#' @return If `len` is scalar, derivative in isotropic lengthscales;
#' if vector, gradient in separable lengthscales.
#' @note require: (XtX, VbtX; k, Jk)
#' @export
gSSDist <- function(len, trunk = 2*k) {
    message(paste(c("[gradient]", "parameter: ", format(len, digits = 8)), collapse = "\t"))
    l <- ifelse(is.matrix(thetaTrain), nrow(thetaTrain), length(thetaTrain))
    np <- length(len)
    K <- kerSEmat(thetaTrain, len = len)
    Kinv <- solve(K)
    lGK <- gradSEmat(thetaTrain, len, K) # not positive semi-definite
    getRGKinv <- function(GK) {
        GKinv <- -Kinv %*% GK %*% Kinv
        GKinv <- forceSymmetric(GKinv)
        rGKinv <- GKinv / Kinv
        rGKinv
    }
    lrGKinv <- lapply(lGK, getRGKinv)
    gd2 <-  vapply(seq_len(l), function(.) gradLOODist2Evd(., Kinv, lrGKinv, VbtX, trunk), double(np))
    if (np == 1) return(sum(gd2))
    return(rowSums(gd2))
}

## VbtX <- with(svdX, tcrossprod(Diagonal(x = d), v)) #0.14s
## lenTest <- seq(0.1, 0.5, by = 0.01)
## listKinv <- map(lenTest, ~getKinv(thetaTrain, .))

## DTscan <- CJ(il = seq(lenTest), ip = seq(thetaTrain))
## DTscan[, len := lenTest[il]]
## DTscan[, theta := thetaTrain[ip]]

## ## SVD version: 11.8s/5 (k=40,l=11)
## system.time(DTscan[, dsvd := map2_dbl(ip, il, ~LOODistSvd(.x, listKinv[[.y]], VbtX))])
## ## Matches raw computation up to 1e-6 relative error.
## ## summary(DTdist$dist / DTscan$dsvd - 1)

## ## EVD version: 37.2s/41 (k=40,l=11); 14.4s/41 (k=20,l=10)
## system.time(DTscan[, devd := map2_dbl(ip, il, ~LOODistEvd(.x, listKinv[[.y]], VbtX))])
## ## Matches raw computation up to 1e-6 relative error, except one with 7.8% rel error.
## ## plot.ecdf(log10(abs(DTdist$dist / DTscan$devd - 1)))

## DTstat2 <- DTscan[, .(e2 = sum(devd^2), e1 = sum(devd)), keyby = .(len)]
## ## lenTest2 <- seq(0.14, 0.25, by = 0.002)
## ## lenTest2 <- setdiff(signif(lenTest2), signif(lenTest))
## ## DTscan2 <- CJ(il = seq(lenTest2), ip = seq(thetaTrain))
## ## DTscan2[, len := lenTest2[il]]
## ## DTscan2[, theta := thetaTrain[ip]]
## ## listKinv <- map(lenTest2, ~getKinv(thetaTrain, .))
## ## DTscan2[, devd := map2_dbl(ip, il, ~LOODistEvd(.x, listKinv[[.y]], VbtX))]
## ## DTstat2 <- rbind(DTstat2, DTscan2[, .(e2 = sum(devd^2), e1 = sum(devd)), keyby = .(len)])
## setkey(DTstat2, len)
## ## k=20,l=10: Both 1-norm and 2-norm approximately selects len=0.19.
## ## k=40,l=11: Both 1-norm and 2-norm approximately selects len=0.25.
## DTstat2[, plot(len, e2, type = 'o', col = "red")]
## DTstat2[, plot(len, e1, type = 'o', col = "black")]
## DTstat2[, plot(len, e2 / max(e2), ylim = c(0.5, 1), type = 'o', col = "red")]
## DTstat2[, lines(len, e1 / max(e1), type = 'o')]
## ## fwrite(DTstat2, paste0("LOOCV-dist-k", k, ".csv"))

## ## Numerical optimization gives len: 0.24833 (k=40); 0.17 (k=20), local mins, e2 smoother than e1.
## ## 14s/17eval; 0.2483433 (default tol = 1.2e-4)
## system.time(ret <- optimize(hSSDist, lower = 0.1, upper = 0.6))
## ## 10.4s/12eval; 0.2470769 (Best option without gradient)
## ##  4.6s/11eval; 0.1682877 (k=20)
## system.time(ret <- optimize(hSSDist, lower = 0.1, upper = 0.6, tol = 1e-2))
## ## 29s/11eval; 0.2483318 (most accurate)
## system.time(ret <- optim(0.2, hSSDist, gSSDist, method = "L-BFGS-B", lower = 0.1, upper = 0.6))
## ## 17s/7eval; 0.2490357 (Best option with gradient)
## ## 10s/8eval; 0.2072725 (k=20; local min!)
## ## 19s/15eval;0.1700722 (k=20; global min; init 0.25)
## tol2 <- list(factr = 10^(15-2))
## tol3 <- list(factr = 10^(15-3))
## system.time(ret <- optim(0.2, hSSDist, gSSDist, method = "L-BFGS-B", lower = 0.1, upper = 0.5,
##                          control = tol3))
