## require(Matrix)

#' Preprocessing via thin SVD
#' @param X Concatenated orthonormal bases.
#' @export
GPSubspacePrepSVD <- function(X) {
    if (!is(X, "Matrix")) X <- Matrix::Matrix(X)
    ## "Big" matrix multiplication: `dpoMatrix`
    XtX <- Matrix::crossprod(X)
    svdX <- svd(X)
    VbtX <- tcrossprod(Matrix::Diagonal(x = svdX$d), svdX$v)
    XtVb <- t(VbtX) # Save the transpose in solving linear equations in the EVD version.
    return(mget(c("XtX", "svdX", "VbtX", "XtVb")))
}
## The combinded basis seems to be full rank (1.129253e-11 > tol = 6.44107e-12)
## r <- Matrix::rankMatrix(X) #0.235s
## message("X has numerical rank ", r, " and ", ncol(X), " columns.")
## system.time(svdX2 <- irlba::irlba(X, t))
## Check that X = u diag(d) t(v)
## with(svdX, all.equal(u %*% tcrossprod(Diagonal(x = d), v), X))

#' Preprocessing via Householder QR with column pivoting
#' @param X Concatenated orthonormal bases.
#' @param transposeR Whether to transpose the R matrix,
#' which saves the transpose in solving linear equations in prediction.
#' @return A list:
#' Pt, transpose of the column pivot matrix.
#' @note R's `base::qr()` uses LINPACK's `DQRDC` for Householder QR with column pivoting.
#' The C++ library Eigen has `ColPivHouseholderQR`, which may be faster
#' and can be accessed via `RcppEigen`.
#' @export
GPSubspacePrepQR <- function(X, transposeR = TRUE) {
    if (!is(X, "Matrix")) X <- Matrix::Matrix(X)
    ## "Big" matrix multiplication: `dpoMatrix`
    XtX <- Matrix::crossprod(X)
    ## QR with column pivoting
    qrX <- qr(X)
    V <- qr.Q(qrX) # `matrix`
    R <- as(qr.R(qrX), "sparseMatrix") # `dgCMatrix`
    Pt <- as(qrX$pivot, "pMatrix") # `pMatrix`
    rankX <- qrX$rank
    if (transposeR) {
        Rt <- t(R)
        return(mget(c("XtX", "V", "Rt", "Pt", "rankX")))
    }
    return(mget(c("XtX", "V", "R", "Pt", "rankX")))
}
