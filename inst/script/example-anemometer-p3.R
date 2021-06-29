## PROM benchmark problem: Anemometer, 3 parameters.
library(gpsr)
library(data.table)
library(purrr)

## Data I/O ----------------------------------------------------------------------
## Objects: Es, Ef, diagonal matrices; Ads, Ac, symmetric; Adf, asymmetric;
## B, a vector, different from 1-parameter anemometer.
## C, a vector, same with 1-parameter anemometer.
## listMat3 <- readRDS("data/anemometer-p3.RDS")
## Data with lazy loading
str(anemometer3)
## Unpack objects into the global environment.
list2env(anemometer3, .GlobalEnv)
n <- nrow(Es)

#' Compute mass matrix
#' @param c Specific heat, [0, 1].
getE <- function(c) {
    Es + c * Ef
}

#' Compute coefficient matrix
#' @param c Specific heat, [0, 1].
#' @param kappa Thermal conductivity, [1, 2].
#' @param v Fluid velocity, [0.1, 2].
getA <- function(c, kappa, v) {
    Ads + kappa * Adf + (c * v) * Ac
}

## Experiment design
d <- 3
nTrain <- 21
nTest <- 100
## Mdesign <- readRDS(paste0("Mdesign", nTrain, ".RDS"))
Mdesign <- Mdesign21

#' Scale a vector in the unit interval to a certain range
#' @param x Numeric vector
#' @param a, b Lower and upper bounds of the new range
toRange <- function(x, a = 0, b = 1) {
    x * (b - a) + a
}
#' Scale a standard design to the anemometer parameter space.
#' @param M Design matrix, 3 columns.
scaleDesign <- function(M) {
    data.table::data.table(c = M[,1],
                           kappa = toRange(M[,2], 1, 2),
                           v = toRange(M[,3], 0.1, 2))
}

DTgrid <- scaleDesign(Mdesign)
pTrain <- DTgrid[seq(nTrain),]
pTest <- DTgrid[-seq(nTrain),]


## Time integration and local POD --------------------------------------------------
timehorizon <- 0.05
timestep <- 1e-3
steps <- as.integer(timehorizon / timestep) #50
k <- 20L

## Compute snapshots
system.time(listE <- purrr::map(DTgrid$c, getE)) #25ms
system.time(listA <- purrr::pmap(DTgrid, getA))  #2s
## system.time(listE <- purrr::map(pTrain$c, getE)) #25ms
## system.time(listA <- purrr::pmap(pTrain, getA))  #2s
system.time(listX <- purrr::map2(listE, listA, ~gpsr:::backwardEuler(.x, .y, B))) #70s/121; 10.7s/21

## Local (raw) POD
system.time(llSVD <- purrr::map(listX, ~irlba::irlba(., k))) #9.4s/118; 1.4s/18
listVpod <- purrr::map(llSVD, ~ .$u)

## Compute reduced system matrices: must provide original parameters without standardization.
system.time(romPOD <- gpsr:::directROM2(DTgrid, listE, listA, listVpod, B, C)) #2.1s; 0.28s
str(romPOD)


## GP subspace regression ------------------------------------------------------------
## (0) Prepare X
listXtrain <- listVpod[seq(nTrain)]
AX <- gpsr:::listMatrix2Array(listXtrain)
str(AX)
l <- nTrain
X <- Matrix::Matrix(AX, n, k * l)

## (1) Overhead computation
system.time(ret <- gpsr::GPSubspacePrepSVD(X)) #2.9s
list2env(ret, .GlobalEnv)

## (2) Hyperparameter tuning
## Require (thetaTrain, XtX, VbtX; Jk)
thetaTrain <- Mdesign[seq(l),]
Jk <- gpsr:::J(k)
## Option 1: Default length-scale, no tuning.
len <- gpsr::defaultLength(d, l)
## Option 2: Use LOOCV error.
lenUpper <- len * 1.3
lenLower <- len * 0.7
system.time(ret <- optimize(gpsr::hSSDist, lower = lenLower, upper = lenUpper, tol = 1e-2)) #4s
## Option 3: Use LOOCV error and gradient.
ToleranceLevel <- function(x) list(factr = 10^(15-x))
tol2 <- ToleranceLevel(2)
system.time(ret2 <- optim(len, gpsr::hSSDist, gpsr::gSSDist, method = "L-BFGS-B",
                          lower = lenLower, upper = lenUpper, control = tol2)) #5.7s
## Not optimal implementation: `optim()` needs separate arguments for objective and gradient,
## while `gSSDist()` can compute both.

## (3) Prediction
len <- ret$minimum
K <- gpsr::kerSEmat(thetaTrain, len = len)
thetaTest <- Mdesign[-seq(l),]
## 2.11s (kl=280); 2.66s (kl=320); 2.91s (kl=340); 3.27s (kl=360); 4.33s/100 (kl=420);
system.time(ret <- purrr::map(seq(nrow(thetaTest)),
                              ~gpsr::GPSubspacePredEVD(thetaTest[.,], thetaTrain, len)))
## (Optional) explicit form of reduced bases: uncessessary in case of (approximate) affine dependency
system.time(listVpred <- purrr::map(ret, ~svdX$u %*% .$Vcirc)) #4.3s

## (4) Compute ROM: must provide original parameters without standardization.
system.time(romGPSR <- gpsr:::directROM2(pTest, listE[-seq(l)], listA[-seq(l)], listVpred, B, C))#1.6s
str(romGPSR)

## Subspace interpolation --------------------------------------------------

## Overhead computation
## Require (listXtrain)
system.time(listGL <- purrr::map(seq(l), gpsr:::GL2i)) #21s;10s;12s;15s;13s
pryr::object_size(listGL) #2GB;0.85GB;1.1GB;1.42GB
## Interpolation
nr <- 5
nTest <- nrow(thetaTest)
system.time(listVint <- purrr::map(seq(nTest),
                                   ~gpsr:::SubspaceInterp(thetaTest[.,], nr, method="RBF")))#10s

## Compute ROM: must provide original parameters without standardization.
system.time(romGInt <- gpsr:::directROM2(pTest, listE[-seq(l)], listA[-seq(l)], listVint, B, C))#1.9s
str(romGInt)


## Error measure: relative L2 state error --------------------------------------------------
#' Relative L2 state error
#' @param i Index of test point
L2StateRelErr <- function(i, listXpred) {
    X <- listX[[i + nTrain]]
    dX <- X - listXpred[[i]] %*% listXr[[i]]
    sqrt(sum(dX^2) / sum(X^2))
}

Ir <- Matrix::Diagonal(x = rep(1, k))

## POD
AAr <- romPOD$AStdAr[,,-seq(l)]
MBr <- romPOD$MStdBr[,-seq(l)]
system.time(listXr <- purrr::map(seq(nTest), ~gpsr:::backwardEuler(Ir, AAr[,,.], MBr[,.])))#3s/100
l2pod <- purrr::map_dbl(seq(nTest), ~L2StateRelErr(., listVpod[-seq(nTrain)]))
summary(l2pod)
mean(l2pod)
##k=21; 18; 14
##8.77e-13; 7.61e-13; 7.98e-13

## GPSR
AAr <- romGPSR$AStdAr
MBr <- romGPSR$MStdBr
system.time(listXr <- purrr::map(seq(nTest), ~gpsr:::backwardEuler(Ir, AAr[,,.], MBr[,.])))#3s/100
l2gpsr <- purrr::map_dbl(seq(nTest), ~L2StateRelErr(., listVpred))
summary(l2gpsr)
mean(l2gpsr)
##8.05e-3; 1.03e-2; 1.47e-2
## +4%, 7.88e-3; +8%, 7.57e-3; +12%, 7.37e-3; +14%, 7.37e-3; +16%, 7.38e-3; +20%, 7.40e-3;

## Subspace interpolation
AAr <- romGInt$AStdAr
MBr <- romGInt$MStdBr
system.time(listXr <- purrr::map(seq(nTest), ~gpsr:::backwardEuler(Ir, AAr[,,.], MBr[,.])))#3s/100
l2gint <- purrr::map_dbl(seq(nTest), ~L2StateRelErr(., listVint))
summary(l2gint)
mean(l2gint)
##2.11e-2; 2.06e-2; 2.55e-2


## Manifold interpolation --------------------------------------------------
## See `interpolate-manifold.R`

## Matrix interpolation --------------------------------------------------
## See `interpolate-matrix.R`
