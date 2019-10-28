fnSampH <- function(s.vec, kappa.2, a.h, b.h.2) {
    s.tilde.vec <- log(s.vec)
    J <- length(s.tilde.vec)
    s.bar <- mean(s.tilde.vec)
    precision <- 1/b.h.2 + J/kappa.2
    draw.mean <- ( (1/b.h.2)*a.h + (J/kappa.2)*s.bar )/precision
    draw.var <- precision^(-1)

    return(rnorm(1, draw.mean, sqrt(draw.var)))
}

fnSampKappa2 <- function(s.vec, h.scal, a.kappa, b.kappa) {
    s.tilde.vec <- log(s.vec)
    J <- length(s.tilde.vec)
    alpha.param <- a.kappa + J/2
    beta.param <- b.kappa + (1/2)*sum( (s.tilde.vec - h.scal)^2)

    return(1/rgamma(1, alpha.param, beta.param))
}
