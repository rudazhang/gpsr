## Numerical schemes to integrate an ODE system.

#' backward Euler (aka implicit Euler) scheme
#' Solves ODE: \eqn{E \dot{x} = A x + b}
#' @param b stationary input vector
#' @return a matrix of state vector snapshots
backwardEuler <- function(E, A, b, dt = timestep, nt = steps, impulse = FALSE) {
    n <- nrow(E)
    X <- matrix(NA_real_, n, nt)
    if (is.matrix(b)) b <- as.vector(b)
    hb <- dt * b
    M <- E - dt * A
    X[,1] <- as.vector(Matrix::solve(M, hb))
    if (!impulse) {
        ## constant input
        for(k in seq(2, nt)) {
            X[,k] <- as.vector(Matrix::solve(M, E %*% X[,k-1] + hb))
        }
    } else {
        ## impulse input
        for(k in seq(2, nt)) {
            X[,k] <- as.vector(Matrix::solve(M, E %*% X[,k-1]))
        }
    }
    return(X)
}

