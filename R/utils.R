## require(Matrix)

## Internal helper functions

#' Vector Euclidean norm
vecnorm <- function(x) sqrt(sum(x^2))

#' Normalize a vector in Euclidean norm
normalize <- function(x) x / sqrt(sum(x^2))

#' Matrix of ones
#' @return A `dpoMatrix` of order n, all entries are 1.
J <- function(n) tcrossprod(Matrix::Matrix(rep(1, n)))

#' Compute the Hadamard product of XtX and kronecker(S, 1_k 1_k^T),
#' without computing the Kronecker product.
#' @note Much slower than XtX * kronecker(DKDinv, Matrix(1, k, k)).
blockHadamard <- function(XtX, S) {
    M <- Matrix::Matrix(NA_real_, k * l, k * l)
    for(i in seq(l)) {
        for(j in seq(l)) {
            rows <- (i-1) * k + seq(k)
            cols <- (j-1) * k + seq(k)
            M[rows, cols] <- XtX[rows, cols] * S[i,j]
        }
    }
    return(M)
}

## #' Whether a `dgCMatrix` is diagonal.
## isDiagonal <- function(M) {
##     ## Row indices of nonzero elements are sequential.
##     identical(M@i, seq(0, nrow(M) - 1)) &
##         ## Each column has one nonzero element, and its a square matrix.
##         identical(M@p, seq(0, nrow(M)))
## }

#' Whether a `dgCMatrix` is a permutation.
isPermutation <- function(M) {
    ## Nonzero elements are all ones.
    all(M@x == 1) &
        ## Each column has one nonzero element, and its a square matrix.
        identical(M@p, seq(0, nrow(M))) &
        ## Row indices of nonzero elements form a permutation
        setequal(M@i, seq(0, nrow(M) - 1))
}

#' Collect a list of matrices into an array.
listMatrix2Array <- function(listMr) {
    nr <- nrow(listMr[[1]])
    nc <- ncol(listMr[[1]])
    vapply(listMr, function(x) as.matrix(x), matrix(NA_real_, nrow = nr, ncol = nc))
}

#' Indicator function of a matrix object
isMatrix <- function(x) {
    is.matrix(x) | is(x, "Matrix")
}

#' Distance vector from a point to a set of points
#' @param x a set of points, either a length-n vector or an n-by-p matrix
#' @param x0 a point, either a scalar or a length-p vector
distvec <- function(x, x0) {
    if (isMatrix(x)) {
        dx <- x - rep(x0, each = nrow(x))
        dist <- apply(dx, 1, vecnorm)
    } else {
        dist <- abs(x - x0)
    }
    return(dist)
}
