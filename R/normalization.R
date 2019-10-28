#' Sample from Dirichlet
#'
#' @export
fnRdirichlet <- function (n, alpha) {
    l <- length(alpha)
    x <- matrix(rgamma(l*n, alpha), ncol=l, byrow=TRUE)
    sm <- x %*% rep(1, l)
    return(x/as.vector(sm))
}

#' Sample Psi
#'
fnSampPsi <- function(d.vec, a.psi) {
    #' @param d.vec vector of length L of cluster sizes
    #' @param a.psi from the prior, vector of length L
    dirichlet.param <- a.psi + d.vec

    return(c(fnRdirichlet(1, dirichlet.param)))
}
