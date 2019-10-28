n.save <- 100
thin <- 1
dat <- readRDS(system.file("extdata", "sim_dat_df.rds", package="bnpzimnr"))

# Settings ----------------------------------------------------------------
SETTINGS <- list()
### MCMC
SETTINGS$print.mcmc <- floor(n.save/10)
### s.j
SETTINGS$s.j.tilde.proposal.sd <- 0.1
### r.ik
SETTINGS$L.r.trunc <- 20
SETTINGS$w.ell.r.proposal.sd <- 1
SETTINGS$r.ik.proposal.sd <- 0.075
### alpha.ij
SETTINGS$L.alpha.trunc <- 150
SETTINGS$w.ell.alpha.proposal.sd <- 1
SETTINGS$alpha.ij.proposal.sd <- 0.15
### beta.kj
SETTINGS$L.beta.trunc <- 50
SETTINGS$beta.kj.proposal.sd <- 0.5
### epsilon.kj
SETTINGS$L.xi.trunc <- 50
SETTINGS$xi.kj.proposal.sd <- 0.3

# Setup -------------------------------------------------------------------
Y.mat <- dat[,-c(1:2)]

r.hat.vec <- log(rowSums(Y.mat)/sum(Y.mat))
r.hat.vec <- r.hat.vec - mean(r.hat.vec)

J <- ncol(Y.mat)
N <- nrow(dat)
n <- length(unique(dat$subject))
K <- N/n
R.rep.mat <- matrix(r.hat.vec, ncol=1)[, rep(1, J)]
Mu.ijk.hat.mat.no.off <- (Y.mat)/exp(R.rep.mat)
y.mat.offset <- 1
Mu.ijk.hat.mat <- (Y.mat+y.mat.offset)/exp(R.rep.mat) # Add offset so the log is not -Inf
Alpha.hat.mat <- log(apply(Mu.ijk.hat.mat, 2, function(otu) tapply(otu, factor(dat$subject), FUN=mean)))
Beta.ikj.hat.mat <- log(Y.mat + y.mat.offset) - R.rep.mat - Alpha.hat.mat
Beta.hat.mat <- apply(Beta.ikj.hat.mat, 2, function(beta.j.vec) tapply(beta.j.vec, factor(dat$condition), FUN=mean))

# Prior -------------------------------------------------------------------
PRIOR <- list()
### s.j
PRIOR$a.s <- 0.25
PRIOR$b.s.2 <- 0.25
### r.ik
PRIOR$a.psi.r <- rep(1.0, SETTINGS$L.r.trunc)
PRIOR$a.w.r <- 5
PRIOR$b.w.r <- 5
PRIOR$u.r.2 <- 0.05
PRIOR$upsilon.r <- mean(r.hat.vec)
PRIOR$b.eta.r.2 <- 0.25
### alpha.ij
PRIOR$a.psi.alpha <- rep(1.0, SETTINGS$L.alpha.trunc)
PRIOR$a.w.alpha <- 1
PRIOR$b.w.alpha <- 1
PRIOR$u.alpha.2 <- 2.0
PRIOR$upsilon.alpha <- mean(log(Mu.ijk.hat.mat.no.off[Y.mat!=0]))
PRIOR$b.eta.alpha.2 <- 1.0
### beta.star.kl
PRIOR$rho.beta <- 1
PRIOR$a.beta.star <- 0
PRIOR$b.beta.star.2 <- 10
PRIOR$a.sigma.beta <- 1
PRIOR$b.sigma.beta <- 1
### xi.star.kl
PRIOR$rho.xi <- 1
PRIOR$a.xi.star <- 2
PRIOR$b.xi.star.2 <- 1
PRIOR$a.sigma.xi <- 1
PRIOR$b.sigma.xi <- 1

# Init --------------------------------------------------------------------
INIT <- list()
## s.j
s.tilde.vec <- rnorm(J, log(0.3), 0.01)
INIT$s.vec <- exp(s.tilde.vec)
## r.ik
INIT$r.vec <- r.hat.vec
INIT$c.r.vec <- .bincode(INIT$r.vec, quantile(INIT$r.vec, (0:5)/5), right=TRUE, include.lowest=TRUE)
INIT$eta.r.vec <- rnorm(SETTINGS$L.r.trunc, 0, 1)
INIT$w.r.vec <- rbeta(SETTINGS$L.r.trunc, 1, 1)
INIT$lambda.r.vec <- rbinom(length(INIT$r.vec), size=1, prob=INIT$w.r.vec[INIT$c.r.vec])
INIT$Delta.mat <- as.integer((Y.mat > 0) + 0)
dim(INIT$Delta.mat) <- c(N, J)
A.mat <- fnGetAMat(INIT$Delta.mat, dat$subject)
B.mat <- fnGetBMat(INIT$Delta.mat, dat$condition)
Delta.sum.mat <- fnGetDeltaSumMat(B.mat)
INIT$Epsilon.mat <- apply(Y.mat, 2, function(otu) tapply(otu==0, factor(dat$condition), FUN=mean))
INIT$Epsilon.mat[INIT$Epsilon.mat<0.01] <- 0.05
INIT$Epsilon.mat[INIT$Epsilon.mat>0.99] <- 0.95
INIT$Epsilon.mat <- INIT$Epsilon.mat + runif(length(INIT$Epsilon.mat), -0.01, 0.01)
dim(INIT$Epsilon.mat) <- c(K, J)
INIT$Xi.mat <- qnorm(INIT$Epsilon.mat)
dim(INIT$Xi.mat) <- dim(INIT$Epsilon.mat)
n.xi.bins <- 5
INIT$C.xi.mat <- t(apply(INIT$Xi.mat, 1, function(xi.k.vec) {
    .bincode(xi.k.vec, quantile(xi.k.vec, (0:n.xi.bins)/n.xi.bins, na.rm=TRUE), right=TRUE, include.lowest=TRUE)
}))
INIT$Xi.star.mat <- matrix(0, nrow=K, ncol=SETTINGS$L.xi.trunc)
for(k in 1:K) {
    xi.star.k.init <- tapply(INIT$Xi.mat[k,], INDEX=INIT$C.xi.mat[k,], FUN=mean)
    INIT$Xi.star.mat[k,] <- replace(INIT$Xi.star.mat[k,], 1:length(xi.star.k.init), xi.star.k.init)
}
INIT$sigma.2.xi.vec <- rep(0.25, K)
# alpha.ij
INIT$Alpha.mat <- Alpha.hat.mat
INIT$c.alpha.vec <- .bincode(INIT$Alpha.mat, quantile(INIT$Alpha.mat, (0:5)/5), right=TRUE, include.lowest=TRUE)
INIT$eta.alpha.vec <- rnorm(SETTINGS$L.alpha.trunc, 0, 1)
INIT$w.alpha.vec <- rbeta(SETTINGS$L.alpha.trunc, 1, 1)
INIT$lambda.alpha.vec <- rbinom(length(INIT$Alpha.mat), size=1, prob=INIT$w.alpha.vec[INIT$c.alpha.vec])
INIT$Alpha.mat[A.mat==0] <- NA
INIT$lambda.alpha.vec[A.mat==0] <- NA
INIT$c.alpha.vec[A.mat==0] <- NA
# beta.kj
INIT$Beta.mat <- Beta.hat.mat
INIT$Beta.mat[Delta.sum.mat==0] <- 0
INIT$Beta.mat[B.mat==0] <- NA
n.beta.bins <- 5
INIT$C.beta.mat <- t(apply(INIT$Beta.mat, 1, function(beta.k.vec) {
    .bincode(beta.k.vec, quantile(beta.k.vec, (0:n.beta.bins)/n.beta.bins, na.rm=TRUE), right=TRUE, include.lowest=TRUE)
}))
INIT$C.beta.mat[INIT$Beta.mat==0] <- NA
INIT$Beta.star.mat <- matrix(0, nrow=K, ncol=SETTINGS$L.beta.trunc)
for(k in 1:K) {
    beta.star.k.init <- tapply(INIT$Beta.mat[k,], INDEX=INIT$C.beta.mat[k,], FUN=mean)
    INIT$Beta.star.mat[k,] <- replace(INIT$Beta.star.mat[k,], 1:length(beta.star.k.init), beta.star.k.init)
}
INIT$sigma.2.beta.vec <- rep(0.25, K)

# Model fit ---------------------------------------------------------------
samples <- zimnr(dat, n.save, thin, INIT, PRIOR, SETTINGS)

