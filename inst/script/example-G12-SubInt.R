## Subspace Interpolation on G_{1,2}.
library(gpsr)
library(data.table)
library(purrr)

thetaTrain <- theta
l <- length(thetaTrain)
k <- 1
listXtrain <- purrr::map(seq(l), ~X[,.])
## thetaTest <- thetanew
thetaTest <- seq(0, 2*pi, length.out = 202)

## Overhead computation: Compute all Grassmann logarithms.
listGL <- purrr::map(seq(l), GL2i)

getG12SubInt <- function(nr) {
    ## Prediction.
    listVint <- purrr::map(thetaTest, ~gpsr:::SubspaceInterp(., nr))
    alpha <- purrr::map_dbl(listVint, ~gpsr:::angleStiefel(.[1,], .[2,]))
    DTsub <- data.table::data.table(theta = thetaTest, alpha = alpha)

    ## Make plot smooth
    DTsub[, t := theta / pi]
    DTsub[, a := alpha / pi]
    DTsub[, delta := c(NA, diff(a))]
    DTsub[, isBranch := delta < 0 & abs((delta + 1) / delta) < 1.5]
    DTsub[1, isBranch := FALSE]
    DTsub[, branch := cumsum(isBranch)]
    DTtail <- DTsub[(isBranch)]
    DTtail[, a := a + 1]
    DTtail[, branch := branch - 1]
    DThead <- DTsub[c(isBranch[-1], FALSE)]
    DThead[, a := a - 1]
    DThead[, branch := branch + 1]
    DTsub <- rbind(DThead, DTsub, DTtail)

    ## Save to file.
    data.table::fwrite(DTsub, paste0("G12-sub-nr", nr, ".csv"))
}

getG12SubInt(nr = 3L)
getG12SubInt(nr = 4L)

