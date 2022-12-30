## Helper functions for the Grassmann manifold G_{1,2} and the MACG distributions on it.

#' Representative angle of a line
#' @param x,y Stiefel representation of a line in plane, vectors of the same length
#' @return a vector of angles, [0, pi)
angleStiefel <- function(x, y) {
    stopifnot(all.equal(x^2 + y^2, rep(1, length(x))))
    alpha <- acos(x * sign(y))
    alpha[y == 0] <- 0
    alpha[]
}

#' PDF of MACG on G_{1,2}, unormalized, maximum as 1.
#' @param a0 Angle relative to the mean; (-pi/2,pi/2].
#' @param c Relative variance, ratio of principal variances: sigma_2^2 / sigma_1^2; [0, 1].
#' @return Probability density of the MACG distribution: p_{MACG}(x; diag(1, c)),
#'         where x = (cos(a0), sin(a0)).
pMACG12 <- function(a0, c) {
    c / (c + (1 - c) * sin(a0)^2)
}

#' Radius of the p-th predictive interval of subspace angle.
#' @param c Relative variance, [0, 1].
#' @param p Probability of predictive interval, (0, 1),
#'          defaults to 0.683 to mimic one standard deviation (68.3%).
alphaPredRadius <- function(c, p = 0.683) {
    atan(sqrt(c) * tan(p / 2 * pi))
}


