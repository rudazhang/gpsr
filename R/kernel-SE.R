#' Squared exponential kernel
#' @param x1,x2 Input parameters, real numbers
#' @param len Correlation length, defaults to pi
## kerSE <- function(x1, x2, len = pi) {
##     exp(-(x1 - x2)^2 / (2 * len^2))
## }

#' Squared exponential kernel, matrix version
#' @param x Vector of scalar parameters, or matrix of vector parameters.
#' @param x0 New parameter. If provided, computes a covariance vector.
#' @param len Length-scale of correlation, isotropic (scalar) or separable (vector), defaults to 1.
#' @return If only x is provided, a `dpoMatrix`
#' (positive-semidefinite matrix from `Matrix`, with Cholesky factor computed).
#' If both are provided, a vector of length nrow(x), if matrix; or length(x), if vector.
#' @export
kerSEmat <- function(x, x0 = NULL, len = 1) {
    ## Isotropic lengthscale
    if(is.matrix(x) & length(len) == 1) {
        len <- rep(len, ncol(x))
    }
    xs <- x %*% Matrix::Diagonal(x = 1 / len)
    if(is.null(x0)) {
        Dist <- as(dist(xs), "Matrix")
        K <- exp(-Dist^2 / 2)
        ## Add a nugget to ensure positive definiteness
        K <- K + Diagonal(ncol(K), 8 * .Machine$double.eps) #4ms
        return(as(K, "dpoMatrix")) # Computes Cholesky factor
    } else {
        x0s <- x0 / len
        dist <- distvec(xs, x0s)
        return(exp(-dist^2 / 2))
    }
}

#' Gradient of the squared exponential kernel, matrix version
#' @param x Vector of scalar parameters, or matrix of vector parameters.
#' @param len Length-scale of correlation, isotropic (scalar) or separable (vector), defaults to 1.
#' @param K Kernel matrix, if already computed
#' @return a list of partial derivatives of the correlation matrix w.r.t. each lengthscale.
gradSEmat <- function(x, len = 1, K = NULL) {
    np <- length(len)
    if(np == 1) {
        Dist <- as(dist(x), "Matrix")
        lFactor <- list(Dist^2 / len^3)
    } else {
        stopifnot(np == ncol(x))
        getFactor <- function(i) {
            Dist <- as(dist(x[,i]), "Matrix")
            Dist^2 / len[i]^3
        }
        lFactor <- lapply(seq(np), getFactor)
    }
    if(is.null(K)) K <- kerSEmat(x, len = len)
    return(lapply(lFactor, function(M) K * M))
}

#' Inverse of correlation matrix
getKinv <- function(x, len) {
    K <- kerSEmat(x, len = len)
    solve(K)
}

#' Rule-of-thumb lengthscale of the SE kernel for GP-subspace
#' @param d Parameter dimension
#' @param l Sample size
#' @export
defaultLength <- function(d, l) 8/3 * sqrt(d) * d / l
