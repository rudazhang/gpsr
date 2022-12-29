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
nTrain <- 18 #21
nTest <- 100
Mdesign <- get(paste0("Mdesign", nTrain))

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
idTrain <- seq(nTrain)
pTrain <- DTgrid[idTrain,]
pTest <- DTgrid[-idTrain,]


## Time integration and local POD --------------------------------------------------
timehorizon <- 0.05
timestep <- 1e-3
steps <- as.integer(timehorizon / timestep) #50
k <- 20L

## Compute snapshots
system.time(listE <- purrr::map(DTgrid$c, getE)) #25ms
system.time(listA <- purrr::pmap(DTgrid, getA))  #2.1s
## system.time(listE <- purrr::map(pTrain$c, getE)) #25ms
## system.time(listA <- purrr::pmap(pTrain, getA))  #2s
system.time(listX <- purrr::map2(listE, listA, ~gpsr:::backwardEuler(.x, .y, B))) #70s/121; 10.7s/21

## Local (raw) POD
system.time(llSVD <- purrr::map(listX, ~irlba::irlba(., k))) #9.4s/118; 1.4s/18
listVpod <- purrr::map(llSVD, ~ .$u)

## Compute reduced system matrices: must provide original parameters without standardization.
system.time(romPOD <- gpsr:::directROM2(DTgrid, listE, listA, listVpod, B, C)) #2.1s; 0.28s
str(romPOD)

## Global POD
AX <- gpsr:::listMatrix2Array(listX[seq(nTrain)])
X <- matrix(AX, n, steps * nTrain)
system.time(retSVD <- irlba::irlba(X, k)) #0.58s; 1.8s
Vpod <- retSVD$u

## Compute reduced system matrices
listVgpod <- purrr::map(seq_along(listE), function(.) return(Vpod))
system.time(romGPOD <- gpsr:::directROM2(DTgrid, listE, listA, listVgpod, B, C)) #2.1s; 0.28s
str(romGPOD)


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
thetaTrain <- Mdesign[idTrain,]
Jk <- gpsr:::J(k)
## Option 1: Default length-scale, no tuning.
(len <- gpsr::defaultLength(d, l))
## Option 2: Use LOOCV error (isotropic lengthscale).
lenUpper <- len * 1.3
lenLower <- len * 0.7
## 4.3s(18); 7.3s (21)
system.time(ret <- optimize(gpsr::hSSDist, lower = lenLower, upper = lenUpper, tol = 1e-2))
## Option 3: Use LOOCV error and gradient (isotropic lengthscale).
ToleranceLevel <- function(x) list(factr = 10^(15-x))
tol2 <- ToleranceLevel(2)
system.time(ret2 <- optim(len, gpsr::hSSDist, gpsr::gSSDist, method = "L-BFGS-B",
                          lower = lenLower, upper = lenUpper, control = tol2)) #4.3s (18); 7s (21)
## Option 4: Use LOOCV error and gradient (separable lengthscale).
## Not an optimal implementation: `optim()` needs separate arguments for objective and gradient,
## while `gSSDist()` can compute both.
## k*l = 420: `hSSDist()`, 1.08s; `gSSDist()`, 1.37s isotropic, 2.03s separable (d=3).
tol3 <- ToleranceLevel(3)
## 12s (18); 23s (21)
system.time(ret3 <- optim(rep(len, d), gpsr::hSSDist, gpsr::gSSDist, method = "L-BFGS-B",
                          lower = rep(0.2, 3), upper = rep(1.5, 3), control = tol3))
## (Option 5: Use LOOCV L2 state error.)

## Summary of lengthscales: l = 21
(len <- gpsr::defaultLength(d, l))
## LOOCV prediction error lengthscales (isotropic gives the best l2 error)
len <- 0.7 ## LOOCV e2: 11.21942; LOOCV l2: 1.726e-2; l2: 8.0e-3;
## 0.008049072 (standardize-then-reduce)
## 0.008003845 (reduce-then-standardize)
len <- c(0.403575, 0.67041, 0.98884)## LOOCV e2: 9.567 (-14.7%); LOOCV l2: 1.85e-2; l2: 9.45e-3 (+18%)
## LOOCV l2 state error lengthscales (isotropic is ok; separable is not good)
len <- 1.6595 ## LOOCV l2: 1.27e-2 (-26.4%); l2: 8.27e-3 (+3.37%);
len <- c(1.635467, 1.407842, 1.797616) ## LOOCV l2: 1.243e-2 (-28.0%); l2: 9.38e-3 (+17%);
## Test l2 state error lengthscales (limited improvement in l2 over LOOCV error isotropic len)
len <- 1.065778 ## LOOCV l2: 1.36e-2 (-21.2%); l2: 7.807e-3 (-2.4%);
len <- c(1.024114, 1.111995, 1.058348) ## LOOCV l2: 1.37e-2 (-20.4%); l2, 7.18e-3 (-10.2%);

## (3) Prediction
len <- ret$minimum
## len <- 0.825 * (1 + 1.5)
K <- gpsr::kerSEmat(thetaTrain, len = len)
thetaTest <- Mdesign[-idTrain,]
## 2.11s (kl=280); 2.66s (kl=320); 2.91s (kl=340); 3.27s (kl=360); 4.33s/100 (kl=420);
system.time(ret <- purrr::map(seq(nrow(thetaTest)),
                              ~gpsr::GPSubspacePredEVD(thetaTest[.,], thetaTrain, len, K)))

## (4) Compute ROM: must provide original parameters without standardization.
## (Optional) explicit form of reduced bases: uncessessary in case of (approximate) affine dependency
system.time(listVpred <- purrr::map(ret, ~svdX$u %*% .$Vcirc)) #4.3s
## Standardize-then-reduce ROMs: 1.8s
system.time(romGPSR <- gpsr:::directROM2(pTest, listE[-idTrain], listA[-idTrain], listVpred, B, C))
## Reduce-then-standardize ROMs; 2.2s
system.time(romGPSR <- gpsr:::directROM2rs(pTest, listE[-idTrain], listA[-idTrain], listVpred, B, C))
str(romGPSR)

## Subspace interpolation --------------------------------------------------

## Overhead computation
## Require (listXtrain)
system.time(listGL <- purrr::map(idTrain, gpsr:::GL2i)) #21s;10s;12s;15s;13s
pryr::object_size(listGL) #2GB;0.85GB;1.1GB;1.42GB
## Interpolation
nr <- 5
nTest <- nrow(thetaTest)
system.time(listVint <- purrr::map(seq(nTest),
                                   ~gpsr:::SubspaceInterp(thetaTest[.,], nr, method="RBF")))#10s

## Compute ROM: must provide original parameters without standardization.
## Standardize-then-reduce ROMs: 1.9s,1.0s
system.time(romGInt <- gpsr:::directROM2(pTest, listE[-idTrain], listA[-idTrain], listVint, B, C))
## Reduce-then-standardize ROMs: 1.4s
system.time(romGInt <- gpsr:::directROM2rs(pTest, listE[-idTrain], listA[-idTrain], listVint, B, C))
str(romGInt)


## Error measure: relative L2 state error --------------------------------------------------
#' Relative L2 state error
#' @param i Index of test point
L2StateRelErr <- function(i, listXpred) {
    X <- listX[[i + nTrain]]
    dX <- X - listXpred[[i]] %*% listXr[[i]]
    sqrt(sum(dX^2) / sum(X^2))
}

Ik <- Matrix::Diagonal(x = rep(1, k))

## POD
AAr <- romPOD$AStdAr[,,-idTrain]
MBr <- romPOD$MStdBr[,-idTrain]
system.time(listXr <- purrr::map(seq(nTest), ~gpsr:::backwardEuler(Ik, AAr[,,.], MBr[,.])))#3s/100
l2pod <- purrr::map_dbl(seq(nTest), ~L2StateRelErr(., listVpod[-seq(nTrain)]))
summary(l2pod)
mean(l2pod)
##k=21; 18; 14
##8.77e-13; 7.61e-13; 7.98e-13

## GPOD
system.time(l2gpod <- TestL2RelErrA(llSVD, romGPOD, listVgpod, timestep, steps))
mean(l2gpod[-seq(nTrain)]) ## l=21a, 5.831e-4 (compare 4.8e-3 of GPS LOOCV L2 isotropic)

## GPSR
AAr <- romGPSR$AStdAr
MBr <- romGPSR$MStdBr
system.time(listXr <- purrr::map(seq(nTest), ~gpsr:::backwardEuler(Ik, AAr[,,.], MBr[,.])))#3s/100
l2gpsr <- purrr::map_dbl(seq(nTest), ~L2StateRelErr(., listVpred))
summary(l2gpsr)
mean(l2gpsr)
## k=21; k=18; k=14
## 8.05e-3; 1.038e-2; 1.465e-2
## k=21: +4%, 7.88e-3; +8%, 7.57e-3; +[12,14]%, 7.37e-3 (best, -8.5%); +16%, 7.38e-3; +20%, 7.40e-3;
## k=18: +40%, 8.45e-3; +80%, 7.57e-3; +120%, 7.23e-3; +150%, 7.17e-3 (best, -31%); +160%, 7.174e-3;
## k=14: -20%, 1.438e-2 (best, -1.8%); -15%, 1.442e-2; -10%, 1.448e-2; -5.7%, 1.455e-2; +5%, 1.47e-2;

## l = 14
## LOOCV error opt (isotropic):
## L2: 0.01464161 (standardize-then-reduce)
## L2: 0.01364237 (reduce-then-standardize)
## LOOCV L2 opt (separable):
## L2: 0.01327473 (standardize-then-reduce)
## L2: 0.01081739 (reduce-then-standardize)

## l = 18
## LOOCV error opt (isotropic):
## len = 0.825
## L2: 0.01037612 (standardize-then-reduce)
## L2: 0.008139351 (reduce-then-standardize)
## len = 2.0625
## L2: 0.007154142 (standardize-then-reduce)
## L2: 0.006486716 (reduce-then-standardize)
## LOOCV L2 opt (separable):


## Subspace interpolation
AAr <- romGInt$AStdAr
MBr <- romGInt$MStdBr
system.time(listXr <- purrr::map(seq(nTest), ~gpsr:::backwardEuler(Ik, AAr[,,.], MBr[,.])))#3s/100
l2gint <- purrr::map_dbl(seq(nTest), ~L2StateRelErr(., listVint))
summary(l2gint)
mean(l2gint)
## k=21; k=18; k=14
## 2.11e-2; 2.06e-2; 2.55e-2

## l = 14
## L2: 0.02552886 (standardize-then-reduce)
## L2: 0.02845684 (reduce-then-standardize)
## 0.01364237 / 0.02845684 = 0.479 (LOOCV error opt, isotropic)
## 0.01081739 / 0.02845684 = 0.380 (LOOCV L2 opt, separable)


## Manifold interpolation --------------------------------------------------
## See `interpolate-manifold.R`

## Matrix interpolation --------------------------------------------------
## See `interpolate-matrix.R`
