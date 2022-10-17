## require(Matrix)
## require(RSpectra)

## Define a small number
## epsMachine <- sqrt(.Machine$double.eps)
## epsMachine <- 1e-12
epsMachine <- .Machine$double.eps

#' GP subspace regression prediction, eigen-decomposition version
#' 36ms, w/o messaging (k=20,l=12); faster than SVD version!
#' @export
GPSubspacePredEVD <- function(thetaNew, thetaTrain, len, K, t = k) {
    l <- ncol(K)
    k <- ncol(XtX) %/% l
    ## message("Predicting at parameter point: ", thetaNew)
    ## Correlation vector
    cv <- kerSEmat(thetaTrain, thetaNew, len)
    ## "Weight" vector: negative component always exists.
    ## DO NOT DO THIS! Using numerical inverse to solve linear equations gives large error.
    ## wv <- Kinv %*% cv
    wv <- solve(K, cv)
    ## DKDinv (aka K_tilde): D_v^-1 K^-1 D_v^-1
    ## Better use linear solve for accuracy. Additional cost is negligible.
    ## wrev <- 1 / wv
    ## DKDinv <- Kinv * tcrossprod(wrev) #1ms
    Dwrev <- Matrix::Diagonal(x = 1 / as.vector(wv))
    Khat <- solve(K, Dwrev)
    DKDinv <- Dwrev %*% Khat #1ms
    ## forceSymmetric() can cause large error. The extra cost of symmpart() is negligible.
    DKDinv <- Matrix::symmpart(DKDinv)
    ## Add a nugget for numerical stability (Must be positive-definite.)
    maxDKDinv <- max(abs(diag(DKDinv)))
    DKDinv <- DKDinv + Matrix::Diagonal(l, epsMachine * maxDKDinv) #4ms
    XPX <- XtX * kronecker(DKDinv, Matrix::Matrix(1, k, k)) #7ms
    ## Add a nugget for numerical stability
    XPX <- XPX + Matrix::Diagonal(k * l, epsMachine * maxDKDinv) #4ms
    ## Force to be a `dpoMatrix`, computes Cholesky factor (upper triangle).
    XPX <- as(XPX, "dpoMatrix") #2ms

    ## Specific to EVD version.
    ## Note: This can be improved by tildeL <- solve(L, XtVb); CPC <- Matrix::crossprod(tildeL),
    ## but `chol()` gives the U factor and t(U) is costly.
    PC <- solve(XPX, XtVb) #3ms
    ## Matrix multiplication: \bar{Sigma} \bar{W}^T U^-1
    CPC <- VbtX %*% PC #6ms
    if (!is.null(t)) { ## Truncated EVD, 2.5ms (almost linear reduction in time)
        evd <- RSpectra::eigs_sym(CPC, t)
    } else { ## Full EVD, 31ms
        evd <- eigen(CPC, symmetric = TRUE)
    }

    ## Noise variance
    eps2 <- 1 - sum(cv * wv)
    if(eps2 < -10 * epsMachine) {
        message("Negative noise variance: ", eps2)
        message("Sample parameters: ", paste(thetaTrain, collapse=", "))
        message("Target parameter: ", thetaNew)
        message("Length scale: ", len)
        stop()
    }
    return(list(Vcirc = evd$vectors, sigma2 = evd$values, eps2 = eps2))
}
## microbenchmark(ret <- GPSubspacePredSVD(0.12, thetaTrain, len, method = "rsvd1")) #49ms; 39ms
## microbenchmark(ret <- GPSubspacePredSVD(0.12, thetaTrain, len, method = "rsvd")) #49ms; 39ms
## microbenchmark(ret <- GPSubspacePredSVD(0.12, thetaTrain, len, method = "svdr")) #43ms; 43ms
## microbenchmark(ret <- GPSubspacePredSVD(0.12, thetaTrain, len, method = "RSpectra")) #51ms; 38ms
## microbenchmark(ret <- GPSubspacePredEVD(0.12, thetaTrain, len)) #35ms; 33ms; 36ms; 36ms; 32ms

