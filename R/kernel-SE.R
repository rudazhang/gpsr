#' Squared exponential kernel
#' @param x1,x2 Input parameters, real numbers
#' @param len Correlation length, defaults to pi
## kerSE <- function(x1, x2, len = pi) {
##     exp(-(x1 - x2)^2 / (2 * len^2))
## }

#' Squared exponential kernel, matrix version
#' @param x Vector of scalar parameters, or matrix of vector parameters.
#' @param x0 New parameter. If provided, computes a covariance vector.
#' @param len Length-scale of correlation, defaults to 1.
#' @return If only x is provided, a `dpoMatrix`
#' (positive-semidefinite matrix from `Matrix`, with Cholesky factor computed).
#' If both are provided, a vector of length nrow(x), if matrix; or length(x), if vector.
#' @export
kerSEmat <- function(x, x0 = NULL, len = 1) {
    if(is.null(x0)) {
        Dist <- as(dist(x), "Matrix")
        K <- exp(-Dist^2 / (2 * len^2))
        ## Add a nugget to ensure positive definiteness
        K <- K + Diagonal(ncol(K), 8 * .Machine$double.eps) #4ms
        return(as(K, "dpoMatrix")) # Computes Cholesky factor
    } else {
        dist <- distvec(x, x0)
        return(exp(-dist^2 / (2 * len^2)))
    }
}

#' Gradient of the squared exponential kernel, matrix version
#' @param K Kernel matrix, if already computed
gradSEmat <- function(x, len = 1, K = NULL) {
    Dist <- as(dist(x), "Matrix")
    Factr <- Dist^2 / len^3
    if(is.null(K)) K <- kerSEmat(x, len = len)
    K * Factr
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
