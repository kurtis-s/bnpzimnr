#' Sample tau^2 vector
#'
#' @export
fnSampTau2Vec <- function(Beta.mat, a.tau, b.tau) {
    J <- ncol(Beta.mat)
    P <- nrow(Beta.mat)

    beta.2.sums <- apply(Beta.mat^2, 1, sum)
    alpha.prime <- a.tau + J/2
    beta.prime <- b.tau + beta.2.sums/2

    tau.2.vec <- 1/rgamma(P, alpha.prime, beta.prime)

    return(tau.2.vec)
}
