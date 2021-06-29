## Utility functions for model reduction.
## require(purrr)

#' Compute reduced standardized coefficient matrix: (E^-1 V)^T A V
reduceStdA <- function(StdV, A, V) Matrix::crossprod(StdV, A %*% V)

#' Compute reduced system matrices given local reduced bases
#' @param p a vector of parameters used for computing the reduced bases
#' @param listA a list of coefficient matrices, from `getA()`.
#' @param listV a list of local reduced bases
#' @return a list: an array of reduced standardized coefficient matrices,
#' matrices of reduce standardized input-to-state vectors and state-to-output vectors.
#' @note Require global objects (for 1-parameter anemometer): E, B, c.
directROM <- function(p, listA, listV, E, B, C) {
    ## "Standardized" reduced bases: E^-1 V
    listStdV <- purrr::map(listV, ~Matrix::solve(E, .)) #0.1s
    ## Reduced standardized coefficient matrix
    listStdAr <- purrr::pmap(list(listStdV, listA, listV), reduceStdA) #0.15s
    AStdAr <- listMatrix2Array(listStdAr) #1ms
    k <- nrow(AStdAr)
    ## Matrix of reduced input-to-state vectors: (E^-1 V)^T B
    MStdBr <- vapply(listStdV, function(.) as.vector(Matrix::crossprod(., B)), double(k)) #0.013s
    ## Matrix of reduced state-to-output vectors: C V
    MCr <- vapply(listV, function(.) as.vector(C %*% .), double(k)) #0.02s
    return(mget(c("p", "AStdAr", "MStdBr", "MCr")))
}

#' @note Require global objects (for 3-parameter anemometer): B, c.
directROM2 <- function(p, listE, listA, listV, B, C) {
    ## "Standardized" reduced bases: E^-1 V
    listStdV <- purrr::map2(listE, listV, ~Matrix::solve(.x, .y)) #0.1s
    ## Reduced standardized coefficient matrix
    listStdAr <- purrr::pmap(list(listStdV, listA, listV), reduceStdA) #0.15s
    AStdAr <- listMatrix2Array(listStdAr) #1ms
    k <- nrow(AStdAr)
    ## Matrix of reduced input-to-state vectors: (E^-1 V)^T B
    MStdBr <- vapply(listStdV, function(.) as.vector(Matrix::crossprod(., B)), double(k)) #0.013s
    ## Matrix of reduced state-to-output vectors: C V
    MCr <- vapply(listV, function(.) as.vector(C %*% .), double(k)) #0.02s
    Mp <- as.matrix(p)
    return(mget(c("Mp", "AStdAr", "MStdBr", "MCr")))
}
