## Document the names of data sets.

#' System matrices for 1-parameter anemometer.
#'
#' System dimension 29,008, with one input and one output.
#' @format A list:
#' \describe{
#'   \item{E}{a diagonal mass matrix}
#'   \item{A1, A2}{sparse coefficient matrices}
#'   \item{B}{input vector}
#'   \item{C}{output (sparse) vector}
#'   \item{P}{a permutation matrix, maps nodes to DOFs, exported from ANSYS}
#' }
#' @source \url{http://modelreduction.org/index.php/Anemometer}
"anemometer1"

#' System matrices for 3-parameter anemometer.
#'
#' System dimension 29,008, with one input and one output.
#' @format A list:
#' \describe{
#'   \item{Es, Ef}{diagonal mass matrices}
#'   \item{Ads, Ac}{symmetric sparse coefficient matrices}
#'   \item{Adf}{asymmetric sparse coefficient matrix}
#'   \item{B}{input vector, different from 1-parameter anemometer}
#'   \item{C}{output (sparse) vector, same with 1-parameter anemometer}
#' }
#' @source \url{http://modelreduction.org/index.php/Anemometer}
"anemometer3"



#' System matrices for thermal model of a micropyros thruster.
#'
#' System dimension 4,257, with one input and seven output; depending on three parameters.
#' @format A list:
#' \describe{
#'   \item{BCIE}{diagonal mass matrix}
#'   \item{BCIA}{symmetric sparse coefficient matrix}
#'   \item{BCIAtop, BCIAbottom, BCIAside}{diagonal coefficient matrices
#'   for the top/bottom/side boundary conditions}
#'   \item{BCIB}{input 1-column sparse matrix}
#'   \item{BCIC}{output 7-row (index) matrix}
#' }
#' @source \url{https://morwiki.mpi-magdeburg.mpg.de/morwiki/index.php/Thermal_Model}
"thermal"
