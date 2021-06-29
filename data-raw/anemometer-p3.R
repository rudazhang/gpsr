## Code to prepare datasets `anemometer-p3.RDS` and `MdesignXX.RDS`.
## PROM benchmark problem: Anemometer, 3 parameters.
## library(lhs)     # Latin hypercube sampling
library(SLHD)    # Maximin-distance Latin hypercube design, a popular LH design.
library(maximin) # batch and sequential maximin design.
library(rgl)
library(FNN)

## Raw Data ----------------------------------------------------------------------
## Objects: E1, E2, diagonal matrices; A1, A2, symmetric; A3-A5, asymmetric;
## B, a column vector, different from 1-parameter anemometer.
## P, a permutation matrix, same with 1-parameter anemometer.
## saveRDS(listMat3, "data/anemometer-p3-raw.RDS")
listMat3 <- readRDS("data/anemometer-p3-raw.RDS")
names(listMat3)

#' Remove small numbers
chop <- function(M, delta = 1e-12) {
    M@x[abs(M@x) < 1e-12] <- 0
    M
}
## Extract the diagonal and form the inverse as a `ddiMatrix`.
Es <- E1
Ef <- E2 - E1
Ads <- A1
Adf <- chop(A3 - A1) + chop(A4 - A5)
Ac <- A2 - A1
## listMat3 <- mget(c("Es", "Ef", "Ads", "Adf", "Ac", "B", "C"))
## saveRDS(listMat3, "param3.RDS")

## Export to MATLAB
## sys <- mget(c("Ads", "Adf", "Ac", "B", "c"))
## sys$Es <- E1
## sys$Ef <- E2 - E1
## system.time(rmatio::write.mat(sys, "anemometer3p.mat", compression = FALSE)) #12ms

## Experiment design ----------------------------------------------------------------------
d <- 3
nTrain <- 18
nTest <- 100

## 1. Generate maximin-LHD sample for training:
## SLHD version: returns Latin square grid; use t = 1 for the standard MmLHD (no slicing).
## lhs version: greedy search, not as good as SLHD.
ret <- SLHD::maximinSLHD(t = 1, m = nTrain, k = d)
str(ret)
## plot(ret$StandDesign, xlim = c(0, 1), ylim = c(0, 1), xaxs = 'i', yaxs = 'i')
rgl::plot3d(ret$StandDesign, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1))

## 2. Generate testing sample:
## 2a. Augment training sample via maximin (preferred): quite slow, but still acceptable.
## Use option for distance to boundary: points clearly away from boundary.
system.time(retAug <- maximin::maximin(nTest, d, Xorig = ret$StandDesign, boundary = TRUE)) #22s
str(retAug)
## points(retAug$Xf[-seq(nTrain),], col = "red")
rgl::open3d()
rgl::plot3d(ret$StandDesign, xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1))
rgl::points3d(retAug$Xf[-seq(nTrain),], col = "red")
## New points are at least as distant from old points than among themselves.
dxaug <- as.vector(FNN::knnx.dist(ret$StandDesign, retAug$Xf[-seq(nTrain),], k = 1))
daug <- as.vector(FNN::knn.dist(retAug$Xf[-seq(nTrain),], 1))
dev.new()
plot.ecdf(dxaug, xlim = c(0, 0.4), xaxs = 'i')
plot.ecdf(daug, add = TRUE, col = "red")

## 2b. Use new maximin-LHD sample: fast; some points are close to the boundary.
retNew <- SLHD::maximinSLHD(t = 1, m = nTest, k = d) #0.95s
str(retNew)
rgl::points3d(retNew$StandDesign, col = "red")
## Half of the new points are closer to old points than among themselves.
dxnew <- as.vector(FNN::knnx.dist(ret$StandDesign, retNew$StandDesign, k = 1))
dnew <- as.vector(knn.dist(retNew$StandDesign, 1))
plot.ecdf(dxnew, xlim = c(0, 0.4), xaxs = 'i')
plot.ecdf(dnew, add = TRUE, col = "red")

## Observation points and target points
Mdesign <- retAug$Xf
saveRDS(Mdesign, paste0("Mdesign", nTrain, ".RDS"))

