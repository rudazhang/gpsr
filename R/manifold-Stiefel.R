## Manifold operations: Stiefel manifold

#' Generate a random sample from the uniform distribution on the Stiefel manifold.
#' @param n Dimension of the Euclidean space.
#' @param k Number of basis vectors.
#' @param seed Seed for the random number generator. Defaults to NA.
rUnifSt <- function(n, k, seed = NA) {
    if (!is.na(seed)) set.seed(seed)
    M <- matrix(rnorm(n * k), n, k)
    return(svd(M)$u)
}

