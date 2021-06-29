## PROM benchmark problem: Anemometer, 1 parameter.
## library(gpsr)
library(Matrix)

## Data I/O ----------------------------------------------------------------------
## listMat1 <- mget(c("E", "A1", "A2", "B", "P", "C"))
## saveRDS(listMat1, "data/anemometer-p1.RDS")

## System matrices:
## E, a diagonal matrix; A1, A2, sparse coefficient matrices; B, c, column vectors.
## P, a permutation matrix, maps nodes to DOFs, exported from ANSYS.
## listMat1 <- readRDS("data/anemometer-p1.RDS")
## Data with lazy loading
str(anemometer1)
## Unpack objects into the global environment.
list2env(anemometer1, .GlobalEnv)
n <- nrow(E)
dA21 <- A2 - A1

#' Compute the coefficient matrix, A
#' @param p Fluid velocity, in [0, 1].
getA <- function(p) A1 + p * dA21

#' Compute the standardized coefficient matrix, E^{-1} A
#' @param p Fluid velocity, in [0, 1].
#' @note  Matrix `%*%` is faster than (A1 + p * dA21) / e.
getStdA <- function(p) solve(E, A1 + p * dA21)

## ROMs by POD ------------------------------------------------------------
k <- 20L
thetaTrain <- round(seq(0, 100, length.out = 7)) / 100
l <- length(thetaTrain)

## Compute snapshots
timehorizon <- 0.05
timestep <- 1e-3    #@Benner2017, Ch. 9
steps <- as.integer(timehorizon / timestep) #50
system.time(listA <- purrr::map(thetaTrain, getA)) #64ms
system.time(listX <- purrr::map(listA, ~backwardEuler(E, ., B, dt = timestep, nt = steps))) #4.2s

## Local POD
system.time(llSVD <- purrr::map(listX, ~irlba::irlba(., k))) #0.55s
listVpod <- purrr::map(llSVD, ~ .$u)

## Compute reduced system matrices
system.time(romPOD <- directROM(thetaTrain, listA, listVpod, E, B, C)) #0.1s
str(romPOD)

## GP subspace regression ------------------------------------------------------------

## (0) Prepare X
listXtrain <- listVpod
AX <- gpsr::listMatrix2Array(listXtrain)
str(AX)
X <- Matrix::Matrix(AX, n, k * l)

## (1) Overhead computation
system.time(ret <- gpsr::GPSubspacePrepSVD(X)) #0.6s
list2env(ret, .GlobalEnv)

## (2) Hyperparameter tuning
## Require (thetaTrain, XtX, VbtX; Jk)
Jk <- J(k)
## Option 1: Default length-scale, no tuning.
len <- gpsr::defaultLength(d = 1, l)
## Option 2: Use LOOCV error.
system.time(ret <- optimize(gpsr::hSSDist, lower = 0.2, upper = 0.8, tol = 1e-2)) #0.5s
## Option 3: Use LOOCV error and gradient.
ToleranceLevel <- function(x) list(factr = 10^(15-x))
tol2 <- ToleranceLevel(2)
system.time(ret2 <- optim(len, gpsr::hSSDist, gpsr::gSSDist, method = "L-BFGS-B",
                          lower = 0.2, upper = 0.8, control = tol2)) #1.2s
## Not optimal implementation: `optim()` needs separate arguments for objective and gradient,
## while `gSSDist()` can compute both.

## (3) Prediction
len <- ret$minimum
## Correlation matrix
K <- gpsr::kerSEmat(thetaTrain, len = len)
thetaTest <- seq(0, 1, by = 0.01)
system.time(ret <- purrr::map(thetaTest, ~gpsr::GPSubspacePredEVD(., thetaTrain, len))) #1s
## (Optional) explicit form of reduced bases: uncessessary in case of (approximate) affine dependency
system.time(listVpred <- purrr::map(ret, ~svdX$u %*% .$Vcirc)) #2s

## (4) Compute ROM
listA <- purrr::map(thetaTest, getA)
romGPSR <- directROM(thetaTest, listA, listVpred, E, B, C)
str(romGPSR)
