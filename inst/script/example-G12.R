## Visualization: GP subspace prediction on G_{1,2}.
library(gpsr)

library(Matrix)
library(data.table)
library(purrr)

## Plot of unormalized PDF of MACG on G_{1,2}.
DTpdf <- data.table::data.table(a0 = seq(-pi/2, pi/2, length.out = 201L))
DTpdf[, p := gpsr:::pMACG12(a0, c = 0.05)]
DTpdf[, plot(a0, p, type = 'l', ylim = c(0, 1), yaxs = 'i')]

## GPS Prediction ----------------------------------------------------------------------
## Observation points and target points
## Seven equi-distance points within 0 and 2 pi.
theta <- seq(0.2, 1.8, length.out = 7) * pi
## Sample size: number of observations
l <- length(theta)

## Test points excluding training points
nNew <- 101
thetanew <- seq(0, 2, length.out = nNew) * pi
thetanew <- setdiff(signif(thetanew), signif(theta))

## Observations: angle (a covering)
alpha <- theta
## Observations: Stiefel representation
X <- rbind(cos(alpha), sin(alpha))
X <- as(X, "dgeMatrix")

## Hyper-parameter: default value of lengthscale.
len <- gpsr::defaultLength(d = 1, l) * (2 * pi)

## GPS Preprocessing
system.time(ret <- gpsr::GPSubspacePrepSVD(X)) #59ms
list2env(ret, .GlobalEnv)
K <- gpsr::kerSEmat(theta, len = len)
## GPS Prediction
predictG12 <- function(target) {
    ret <- gpsr::GPSubspacePredEVD(target, theta, len, K, t = NULL)
    v <- svdX$u %*% ret$Vcirc[,1]
    alpha <- gpsr:::angleStiefel(v[1,], v[2,])
    lambda <- ret$sigma2 + ret$eps2
    c <- lambda[[2]] / lambda[[1]]
    c(alpha = alpha, c = c)
}
system.time(Mpred <- vapply(thetanew, predictG12, double(2)))
DTpred <- data.table::as.data.table(t(Mpred))
DTpred$radius <- gpsr:::alphaPredRadius(DTpred$c, p = 0.95)
DTpred$theta <- thetanew
DTtrain <- data.table::data.table(alpha = theta %% pi, c = 0, radius = 0, theta = theta)
DTpred <- rbind(DTpred, DTtrain)
data.table::setkey(DTpred, theta)

## Basic plot
DTpred[, plot(theta, alpha, type = 'l', ylim = c(0, pi), yaxs = 'i')]
points(theta, theta %% pi, pch = 19)
curve(x %% pi, col = "blue", add = TRUE)
DTpred[, lines(theta, alpha + radius, col = "red")]
DTpred[, lines(theta, alpha - radius, col = "red")]


## Plots --------------------------------------------------------------------------------
## Make plot smooth (and normalized)
DTpred[, c("t", "a", "r") := lapply(.SD, `/`, pi), .SDcols = c("theta", "alpha", "radius")]
DTpred[, delta := c(0, diff(a))]
DTpred[, dist := delta %% 1]
DTpred[, branch := dist - delta]
DTpred[, correction := cumsum(branch)]
DTpred[, a0 := a + correction]

## Subspace interpolation results.
if (!file.exists("G12-sub-nr3.csv")) {
    source("inst/script/example-G12-SubInt.R")
}
DTsub3 <- data.table::fread("G12-sub-nr3.csv")
DTsub4 <- data.table::fread("G12-sub-nr4.csv")

## 2D plot
colTruth <- "black"
colGPSR <- "blue"
colSubspace <- "orange"
shadeColor <- adjustcolor("red", 0.2)
textCex <- 2
lineWidth <- 2

png("G12-plot.png", width = 1000, height = 600)
## png("G12-plot-small.png", width = 700, height = 350)
## lineWidth <- 1.5
op <- par(mar = c(4, 5, 1, 1), las = 1)
DTpred[, plot(t, a, type = 'n', xlim = c(0, 2), xaxs = 'i', ylim = c(0, 1), yaxs = 'i',
              xlab = "", ylab = "", yaxt = 'n', cex.axis = textCex)]
axis(2, at = seq(0, 1, by = 0.5), tick = FALSE, cex.axis = textCex)
axis(2, at = seq(0, 1, by = 0.25), labels = FALSE, cex.axis = textCex)
mtext("parameter (pi)", 1, line = 2.5, cex = textCex)
par(las = 0)
mtext("subspace angle (pi)", 2, line = 3.5, cex = textCex)
## Subspace interpolation
DTsub3[, lines(t, a, col = colSubspace, lwd = lineWidth), keyby = .(branch)]
DTsub4[, lines(t, a, col = colSubspace, lwd = lineWidth, lty = 3), keyby = .(branch)]
## True function
abline(a = 0, b = 1, col = colTruth, lwd = lineWidth)
abline(a = -1, b = 1, col = colTruth, lwd = lineWidth)
## Sample points
points(theta / pi, (theta / pi) %% 1, pch = 19)
#' Add to the plot one continuous branch of angle representation of the probabilistic prediction.
#' @param branch Integer of branch offset.
AddBranchPrediction <- function(off = 0) {
    DTplot <- DTpred[, .(t, a = a0 + off, r)]
    DTplot[, lines(t, a, col = colGPSR, lwd = lineWidth)]
    DTplot[, polygon(c(t, rev(t)), c(a - r, rev(a + r)), col = shadeColor, border = NA)]
}
AddBranchPrediction(0)
AddBranchPrediction(-1)
AddBranchPrediction(-2)
AddBranchPrediction(1)
## legend(0.5, 0.5, bty = "n", lwd = 2, cex = textCex,
##        legend = c("true map", "GP-subspace", "Subspace-Int, nr=3", "Subspace-Int, nr=4"),
##        col = c("black", colGPSR, colSubspace, colSubspace),
##        lty = c(1, 1, 1, 3))
par(op)
dev.off()

## 2D Plot: raster
pointCex <- 2
lineWidth <- 2

png("G12-raster.png", width = 1400, height = 700)
op <- par(mar = c(0,0,0,0))
DTpred[, plot(t, a, type = 'n', xlim = c(0, 2), xaxs = 'i', ylim = c(0, 1), yaxs = 'i',
              bty = "n", xaxt = "n", yaxt = "n",
              xlab = "theta / pi", ylab = "alpha / pi")]
## True function
abline(a = 0, b = 1, col = colTruth, lwd = lineWidth)
abline(a = -1, b = 1, col = colTruth, lwd = lineWidth)
## Subspace interpolation
DTsub3[, lines(t, a, col = colSubspace, lwd = lineWidth), keyby = .(branch)]
DTsub4[, lines(t, a, col = colSubspace, lwd = lineWidth, lty = 3), keyby = .(branch)]
## GPS prediction
AddBranchPrediction(0)
AddBranchPrediction(-1)
AddBranchPrediction(-2)
AddBranchPrediction(1)
## Sample points
points(theta / pi, (theta / pi) %% 1, pch = 19, cex = pointCex)
par(op)
dev.off()

## Double the plot with imagemagick.
hasConvert <- (system('which convert', ignore.stdout = TRUE) == 0L)
if (hasConvert) {
    system('convert -append G12-raster.png G12-raster.png G12-raster-2pi.png')
} else {
    library(magick)
    imgpi <- magick::image_read("G12-raster.png")
    img2pi <- magick::image_append(c(imgpi, imgpi), stack = TRUE)
    magick::image_write(img2pi, "G12-raster-2pi.png")
}

## 3D plot: raster 2d plot as texture on 3d cylinder
library(rgl)
## Define a cylinder.
sides <- 40L
nodes <- sides + 1L
Theta <- matrix(c(0, 2 * pi), 2, nodes)
Alpha <- matrix(seq(0, 2 * pi, length.out = nodes), 2, nodes, byrow = TRUE)
x <- Theta
y <- cos(Alpha)
z <- sin(Alpha)

## Make 3d plot.
persp3d(x, y, z, aspect = "iso",
        texture = "G12-raster-2pi.png",
        specular = "white", col = "white",
        axes = FALSE, box = FALSE, xlab = "", ylab = "", zlab = "")
lines3d(x = c(0, 2*pi), y = 0, z = 0, col = "black", lwd = 3)
## Saved viewpoint
par3d(userMatrix = matrix(c(0.89381659, 0.3121026, 0.3219992, 0,
                            -0.03827326, -0.6623411, 0.7482244, 0,
                            0.44679594, -0.6810993, -0.5800664, 0,
                            0, 0, 0, 1), nrow = 4, byrow = TRUE),
      windowRect = c(10, 82, 970, 1125), zoom = 1)
## SVG export has no texture image.
## rgl.postscript("G12-rgl.svg", fmt = "svg")
snapshot3d("G12-rgl.png", webshot = FALSE)

## Trim white margins.
if (hasConvert) {
    system('convert G12-rgl.png -trim G12-rgl-trim.png')
} else {
    img <- magick::image_read("G12-rgl.png")
    imgtrim <- magick::image_trim(img)
    magick::image_write(imgtrim, "G12-rgl-trim.png")
}
