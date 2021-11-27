## Manifold operations
## require(Matrix)

#' Grassmann logarithm
#' @param X0 Stiefel representation of base point on the Grassmann manifold
#' @param Y Stiefel representation of target point
#' @return the horizontal lift of log_[X0] [Y]
#' @references @Bendokat2020
GrassmannLog <- function(X0, Y) {
    ## Y^T X = Q S R^T
    svdTilde <- svd(Matrix::crossprod(Y, X0))
    ## bar{Y} = Y (Q R^T)
    barY <- Y %*% tcrossprod(svdTilde$u, svdTilde$v)
    ## bar{X} = X (R S R^T)
    barX <- X0 %*% (svdTilde$v %*% Matrix::tcrossprod(Matrix::Diagonal(x = svdTilde$d), svdTilde$v))
    ## bar{Y} - bar{X} = Q S R^T
    svdHat <- svd(barY - barX)
    ## Delta_x = Q arcsin(S) R^T
    Delta <- with(svdHat, u %*% Matrix::tcrossprod(Matrix::Diagonal(x = asin(d)), v))
    Delta
}
## GL <- GrassmannLog(X0, Y)
## Classical method to compute the Grassmann logarithm.
## box0i <- solve(Matrix::crossprod(X0, Y))
## ret <- svd(Y %*% box0i - X0)
## GL2 <- with(ret, u %*% tcrossprod(Diagonal(x = atan(d)), v))
## all.equal(GL, GL2)

#' Grassmann exponential
#' @param T horizontal lift of a tangent vector at X0.
#' @return a Stiefel representation of exp_[X0](T)
GrassmannExp <- function(X0, T) {
    svdT <- svd(T)
    ret <- with(svdT, X0 %*% (v %*% Matrix::Diagonal(x = cos(d))) +
                      u %*% Matrix::Diagonal(x = sin(d)))
    ret
}

#' Principal angles between two subspaces
#' @return a vector of principal angles in increasing order, [0, pi/2]
PrAngles <- function(X, Y) {
    XtY <- Matrix::crossprod(X, Y)
    sigma <- svd(XtY, nu = 0, nv = 0)$d
    ## suppress numerical errors for singular values near 1.
    sigma[sigma > 1] <- 1
    acos(sigma)
}

#' The Riemannian distance between two subspaces: 2-norm of principal angles
RiemannianDist <- function(X, Y, normalize = TRUE) {
    angles <- PrAngles(X, Y)
    if (normalize) angles <- angles / (pi / 2)
    sqrt(sum(angles^2))
}

#' The gap metric between two subspaces: sine of the largest principal angle
#' @param X, Y Stiefel representations.
GapMetric <- function(X, Y) {
    norm(Matrix::crossprod(X) - Matrix::crossprod(Y))
}

#' Angle between a vector and a subspace
#' @param X an orthonormal frame of a subspace, an n-by-k matrix
#' @return a scalar, [0, pi/2]
AngleVS <- function(v, X) {
    coef <- Matrix::crossprod(X, v)
    cosAngle <- vecnorm(coef) / vecnorm(v)
    ifelse(cosAngle > 1, 0, acos(cosAngle))
}
