## Utility functions for model reduction.
## require(purrr)

#' Galerkin projection
galerkin <- function(M, V) Matrix::crossprod(V, M %*% V)

#' Compute reduced standardized coefficient matrix: (E^-1 V)^T A V
#' @note standardize-then-reduce
reduceStdA <- function(StdV, A, V) Matrix::crossprod(StdV, A %*% V)

#' Compute reduced system matrices given local reduced bases.
#'
#' Single-input single-output, varying coefficient matrix, Galerkin projection.
#' @param p a vector of parameters used for computing the reduced bases
#' @param listA a list of coefficient matrices, from `getA()`.
#' @param listV a list of local reduced bases
#' @return a list: an array of reduced standardized coefficient matrices,
#' matrices of reduce standardized input-to-state vectors and state-to-output vectors.
#' @note Require global objects (for 1-parameter anemometer): E, B, c.
#' @note standardize-then-reduce
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

#' Varying mass and coefficient matrices
#' @note standardize-then-reduce
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

#' Compute reduced system matrices given local reduced bases.
#' @note reduce-then-standardize
directROMrs <- function(p, listA, listV, E, B, C) {
    ## Reduced system matrices
    listEr <- purrr::map(listV, ~galerkin(E, .))
    listAr <- purrr::map2(listA, listV, galerkin)
    listBr <- purrr::map(listV, ~Matrix::crossprod(., B))
    ## Standardized reduced system matrices
    listStdAr <- purrr::map2(listEr, listAr, ~Matrix::solve(.x, .y))
    listStdBr <- purrr::map2(listEr, listBr, ~Matrix::solve(.x, .y))
    AStdAr <- listMatrix2Array(listStdAr) #1ms
    k <- nrow(AStdAr)
    ## Matrix of reduced input-to-state vectors
    MStdBr <- vapply(listStdBr, function(.) as.vector(.), double(k))
    ## Matrix of reduced state-to-output vectors: C V
    MCr <- vapply(listV, function(.) as.vector(C %*% .), double(k)) #0.02s
    Mp <- as.matrix(p)
    return(mget(c("Mp", "AStdAr", "MStdBr", "MCr")))
}

#' Varying mass and coefficient matrices
#' @note reduce-then-standardize
directROM2rs <- function(p, listE, listA, listV, B, C) {
    ## Reduced system matrices
    listEr <- purrr::map2(listE, listV, galerkin)
    listAr <- purrr::map2(listA, listV, galerkin)
    listBr <- purrr::map(listV, ~Matrix::crossprod(., B))
    ## Standardized reduced system matrices
    listStdAr <- purrr::map2(listEr, listAr, ~Matrix::solve(.x, .y))
    listStdBr <- purrr::map2(listEr, listBr, ~Matrix::solve(.x, .y))
    AStdAr <- listMatrix2Array(listStdAr) #1ms
    k <- nrow(AStdAr)

    ## Matrix of reduced input-to-state vectors
    MStdBr <- vapply(listStdBr, function(.) as.vector(.), double(k))
    ## Matrix of reduced state-to-output vectors: C V
    MCr <- vapply(listV, function(.) as.vector(C %*% .), double(k)) #0.02s
    Mp <- as.matrix(p)
    return(mget(c("Mp", "AStdAr", "MStdBr", "MCr")))
}

#' Compute reduced system matrices given a global reduced basis.
#' @note standardize-then-reduce
directGROM <- function(p, listA, V, E, B, C) {
    ## "Standardized" reduced bases: E^-1 V
    StdV <- Matrix::solve(E, V) #0.1s
    ## Reduced standardized coefficient matrix
    listStdAr <- purrr::map(listA, ~reduceStdA(StdV, ., V)) #0.15s
    AStdAr <- listMatrix2Array(listStdAr) #1ms
    ## Matrix of reduced input-to-state vectors: (E^-1 V)^T B
    l <- length(listA)
    StdBr <- as.vector(Matrix::crossprod(StdV, B))
    MStdBr <- matrix(rep(StdBr, l), ncol = l)
    ## Matrix of reduced state-to-output vectors: C V
    Cr <- as.vector(C %*% V)
    MCr <- matrix(rep(Cr, l), ncol = l)
    return(mget(c("p", "AStdAr", "MStdBr", "MCr")))
}

#' Compute reduced system matrices given a global reduced basis.
#' @note reduce-then-standardize
directGROMrs <- function(p, listA, V, E, B, C) {
    ## Reduced system matrices
    Er <- galerkin(E, V)
    listAr <- purrr::map(listA, ~galerkin(., V))
    Br <- Matrix::crossprod(V, B)
    ## Standardized reduced system matrices
    listStdAr <- purrr::map(listAr, ~Matrix::solve(Er, .))
    AStdAr <- listMatrix2Array(listStdAr) #1ms
    ## Matrix of reduced input-to-state vectors
    l <- length(listA)
    StdBr <- as.vector(Matrix::solve(Er, Br))
    MStdBr <- matrix(rep(StdBr, l), ncol = l)
    ## Matrix of reduced state-to-output vectors: C V
    Cr <- as.vector(C %*% V)
    MCr <- matrix(rep(Cr, l), ncol = l)
    return(mget(c("p", "AStdAr", "MStdBr", "MCr")))
}
