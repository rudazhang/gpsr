## require(Matrix)
## require(RSpectra)
## require(rsvd)

## NOTE: Remove `system.time()` and messaging when not debugging a function!
#' 44ms, w/o messaging (k=20,l=12).
#' GP subspace regression prediction, posterior from equal likelihood
#' @note Require (K, XtX, VbtX)
#' @param thetaNew a target parameter point
#' @param t        truncation level, defaults to k; if NULL, use thin SVD.
#' @return a list: Vcirc, principal directions as coordinates in the global basis,
#'                 a c-by-t matrix (more than the mean prediction);
#'                 sigma2, principal variances, a vector of length t;
#'                 eps2, noise variance, a scalar.
#' @details To get the explicit principal directions, do matrix multiplication V = stdX$u %*% Vcirc.
#' @export
GPSubspacePredSVD <- function(thetaNew, thetaTrain, len, t = k, method = "rsvd") {
    l <- ncol(K)
    k <- ncol(XtX) %/% l
    ## message("Predicting at parameter point: ", thetaNew)
    ## Correlation vector
    cv <- kerSEmat(thetaTrain, thetaNew, len)
    ## "Weight" vector: negative component always exists.
    ## DO NOT DO THIS! Using numerical inverse to solve linear equations gives large error.
    ## wv <- Kinv %*% cv
    wv <- solve(K, cv)
    ## DKDinv (aka Pi): D_v^-1 K^-1 D_v^-1
    ## Better use linear solve for accuracy. Additional cost is negligible.
    ## wrev <- 1 / wv
    ## DKDinv <- Kinv * tcrossprod(wrev) #1ms
    Dwrev <- Matrix::Diagonal(x = 1 / as.vector(wv))
    Khat <- solve(K, Dwrev)
    DKDinv <- Dwrev %*% Khat #1ms
    ## forceSymmetric() can cause large error. The extra cost of symmpart() is negligible.
    DKDinv <- Matrix::symmpart(DKDinv)
    ## Add a nugget for numerical stability
    maxDKDinv <- max(abs(diag(DKDinv)))
    DKDinv <- DKDinv + Matrix::Diagonal(l, epsMachine * maxDKDinv) #4ms
    XPX <- XtX * kronecker(DKDinv, Matrix::Matrix(1, k, k)) #7ms
    ## Add a nugget for numerical stability
    XPX <- XPX + Matrix::Diagonal(k * l, epsMachine * maxDKDinv) #4ms
    ## Force to be a `dpoMatrix`, computes Cholesky factor (upper triangle).
    XPX <- as(XPX, "dpoMatrix") #2ms

    ## Specific to SVD version.
    U <- chol(XPX) #0ms
    Uinv <- solve(U) #3ms
    ## Matrix multiplication: \bar{Sigma} \bar{W}^T U^-1
    SWUinv <- VbtX %*% Uinv #6ms
    if (!is.null(t)) { ## Truncated SVD: 6ms (k=20)
        svdSmall <- rsvd::rsvd(SWUinv, t)
    } else { ## Thin SVD: 50ms (r=240)
        svdSmall <- svd(SWUinv)
    }

    ## Noise variance
    eps2 <- 1 - sum(cv * wv)
    if(eps2 < 0) {
        message("Negative noise variance: ", eps2)
        message("Sample parameters: ", paste(thetaTrain, collapse=", "))
        message("Target parameter: ", thetaNew)
        message("Length scale: ", len)
        stop()
    }
    return(list(Vcirc = svdSmall$u, sigma2 = svdSmall$d^2, eps2 = eps2))
}
## Truncated SVD vs. truncated EVD
## system.time(svdSmall <- svd(SWUinv)) #702ms; 3.03s
## microbenchmark(svdSmall1 <- irlba::irlba(SWUinv, t), times = 100)   #21ms; 370ms
## microbenchmark(svdSmall2 <- irlba::svdr(SWUinv, t), times = 100)    #10ms; 159ms; 914ms
## microbenchmark(svdSmall3 <- RSpectra::svds(SWUinv, t), times = 100) #8ms; 31.3ms; 81.4ms
## microbenchmark(svdSmall4 <- rsvd::rsvd(SWUinv, t), times = 100)     #6ms; 104ms (k*l=840); 193ms
## microbenchmark(svdSmall5 <- rsvd::rsvd(SWUinv, t, q = 1), times = 100) #4ms; 69ms (k*l=840); 156ms
## system.time(evd <- eigen(CPC, symmetric = TRUE)) #417ms; 1.23s
## microbenchmark(evd <- RSpectra::eigs_sym(CPC, t), times = 100)      #3ms; 12.6ms (k*l=840); 22.2ms
