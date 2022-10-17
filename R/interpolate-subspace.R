## ROMs by subspace interpolation
## NOTE: Principal angles mostly grow linearly with parameter,
## which verifies the premise of subspace interpolation!

#' Compute Grassmann logarithm at p_i for all points p_j, j != il.
#' @note Require (listXtrain)
GL2i <- function(i) {
    l <- length(listXtrain)
    id <- setdiff(seq(l), i)
    ret <- purrr::map(id, ~GrassmannLog(listXtrain[[i]], listXtrain[[.]]))
    names(ret) <- id
    ret
}

#' Tangent vector Lagrange interpolation (one-parameter)
#' @param theta New parameter point, a scalar
#' @param ref   Index of the reference training parameter point
#' @param neigh Indices of the used neighboring training parameter points
#' @return Tangent vector
#' @note Require (thetaTrain, listGL)
TangentInterpLagrange <- function(theta, ref, neigh) {
    nn <- length(neigh)
    ## Weight vector for Lagrange interpolation
    LagrangeWeight <- LagrangeInterp(thetaTrain[c(ref, neigh)])
    weights <- LagrangeWeight(theta)
    cNeigh <- as.character(neigh)
    i <- 1L
    ret <- listGL[[ref]][[cNeigh[i]]] * weights[i+1]
    for(i in seq(2, nn)) {
        ret <- ret + listGL[[ref]][[cNeigh[i]]] * weights[i+1]
    }
    return(ret)
}
## ref <- 2
## neigh <- c(1, 3, 4)
## theta <- 0.12

#' Tangent vector interpolation with (multiquadric) radial basis function
#' @param theta New parameter point, a length-d vector or a scalar
#' @note Require (thetaTrain, listGL, listXtrain)
TangentInterpRBF <- function(theta, ref, neigh) {
    X0 <- listXtrain[[ref]]
    n <- nrow(X0)
    k <- ncol(X0)
    nn <- length(neigh)
    cNeigh <- as.character(neigh)
    zerotan <- double(n * k)
    Mtan <- vapply(seq(nn), function(i) as.vector(listGL[[ref]][[cNeigh[i]]]), zerotan)
    Mtan <- cbind(matrix(zerotan, ncol = 1), Mtan)
    if (is.matrix(thetaTrain)) {
        pInt <- thetaTrain[c(ref, neigh), ]
    } else {
        pInt <- thetaTrain[c(ref, neigh)]
    }
    Interp <- MultiquadricRBF(pInt, t(Mtan))
    M <- matrix(Interp(theta), nrow = n)
    ## Horizontal projection
    Z <- M - X0 %*% crossprod(X0, M)
    return(Z)
}

#' Subspace interpolation
#' @param theta New parameter point
#' @param nr Number of (parameter, reduced basis) pairs used in the interpolation.
#'           Defaults to all observations.
#' @param method Interpolation scheme, either "Lagrange" or "RBF" (for multiquadrics)
#' @return An orthonormal basis for the interpolated subspace.
#' @references @Amsallem2008
#' @note Require (thetaTrain, listXtrain)
#' @details It not clear from [@Amsallem2008] and later papers that
#' whether the reference point should be one of the interpolation points,
#' but it is more reasonable to do than exclude it.
#' The tangent vector of the reference point is identically zero.
SubspaceInterp <- function(theta, nr = l, method = "Lagrange") {
    ## Need at least 3 nearby points.
    ## If nr = 2, use piecewise linear interpolation (not implemented).
    stopifnot(nr >= 3)
    dist <- gpsr:::distvec(thetaTrain, theta)
    id <- order(dist)[seq(nr)]
    ref <- id[1]
    neigh <- id[-1]
    if (method == "Lagrange") {
        Tangent <- TangentInterpLagrange(theta, ref, neigh)
    } else if (method == "RBF") {
        Tangent <- TangentInterpRBF(theta, ref, neigh)
    } else {
        stop("`method` must be 'Lagrange' or 'RBF'. Other schemes are not implemented.")
    }
    GrassmannExp(listXtrain[[ref]], Tangent)
}
