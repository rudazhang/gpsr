## Interpolation Schemes

#' A function factory for Lagrange interpolation.
#' @param p Interpolation points
#' @return a function that takes in a new point and
#' returns the weight vector for Lagrange interpolation.
LagrangeInterp <- function(p) {
    l <- length(p)
    P <- tcrossprod(p, rep(1, l))
    DP <- P - t(P)
    diag(DP) <- 1
    pi <- apply(DP, 1, prod)
    rm(P, DP)
    function(x) {
        dx <- x - p
        pix <- vapply(seq_along(p), function(i) prod(dx[-i]), double(1))
        pix / pi
    }
}
## p <- seq(0, 1, length.out = 2)
## np <- length(p)
## LagrangeWeight <- LagrangeInterp(p)
## x <- seq(-.3, 1.3, by = 0.01)
## Mw <- vapply(x, LagrangeWeight, double(np))
## str(Mw)
## matplot(x, t(Mw), type = 'l', lty = 1)
## abline(h = c(0, 1), lty = 2)
## points(p, rep(1, np), pch = 19, col = seq(np))

#' A function factory for piecewise linear interpolation.
#' @param p Interpolation points, can be unordered.
#' @return a function that takes in a new point and
#' returns the weight vector for Lagrange interpolation.
PiecewiseLinearInterp <- function(p) {
    l <- length(p)
    ord <- order(p)
    ps <- sort(p)
    function(x) {
        w <- double(l)
        idLess <- which(x > ps)
        if (length(idLess) == 0) {
            w[ord[1]] <- 1
            return(w)
        }
        a <- max(idLess)
        if (a == l) {
            w[ord[l]] <- 1
            return(w)
        }
        b <- a + 1
        pa <- ps[a]
        pb <- ps[b]
        wa <- (pb - x) / (pb - pa)
        wb <- 1 - wa
        w[ord[a]] <- wa
        w[ord[b]] <- wb
        w
    }
}
## p <- seq(0,1, by = 0.25)
## np <- length(p)
## PiecewiseLinearWeight <- PiecewiseLinearInterp(p)
## x <- seq(0, 1.1, by = 0.01)
## Mw <- vapply(x, PiecewiseLinearWeight, double(np))
## str(Mw)
## matplot(x, t(Mw), type = 'l', lty = 1)
## abline(h = c(0, 1), lty = 2)
## points(p, rep(1, np), pch = 19, col = seq(np))

#' Hardy's multiquadrics, a radial basis function
#' @param p Interpolation points, a n-by-d matrix
#' @param f Function values, a length-n vector or an n-by-q matrix (vector-valued function)
#' @return A function that takes in a new point
#' and returns the value of the interpolation function.
#' @note Parameter computed per [@Amsallem2010, Appendix H]
MultiquadricRBF <- function(p, f) {
    ## d <- ncol(p)
    ranges <- apply(p, 2, function(x) diff(range(x)))
    maxRange <- max(ranges)
    R <- sqrt(maxRange / 10)
    q <- 0.25
    D <- as.matrix(dist(p))
    D <- as(D, "dsyMatrix")
    B <- (D^2 + R^2)^q
    a <- solve(B, f)
    function(x) {
        dist <- distvec(p, x)
        b <- (dist^2 + R^2)^q
        as.vector(crossprod(b, a))
    }
}
## ## 1d test case: sine curve
## p <- seq(0, 2 * pi, length.out = 7)
## f <- sin(p)
## H <- MultiquadricRBF(matrix(p, ncol = 1), f)
## x <- seq(0, 2 * pi, length.out = 120)
## h <- vapply(x, H, double(1))
## ## Compare with Lagrange interpolation
## HL <- LagrangeInterp(p)
## wL <- vapply(x, HL, double(length(p)))
## hL <- crossprod(wL, f)
## plot(x, h, type = 'l')
## points(p, f)
## lines(x, hL, col = "red")
## ## vector-valued test case: sine and cosine curves
## f <- cbind(sin(p), cos(p))
## H <- MultiquadricRBF(matrix(p, ncol = 1), f)
## h <- vapply(x, H, double(2))
## plot(x, h[1,], type = 'l', col = "red")
## points(p, f[,1], col = "red")
## lines(x, h[2,], col = "blue")
## points(p, f[,2], col = "blue")
