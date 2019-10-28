#' Sample stick-breaking weights
#'
#' @param rho DP concentration parameter
#' @export
fnSampVVec <- function(rho, d.vec, d.cumsum.vec) {
    L <- length(d.vec)
    rbeta(L-1, 1+d.vec[-L], rho + d.cumsum.vec)
}
