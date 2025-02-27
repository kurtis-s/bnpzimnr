#' Matrix with each element a_{ij} = I(\\sum_{k=1}^K \\delta_{ikj} \\ne 0)
#'
#' @export
fnGetAMat <- function(Delta.mat, subject) {
    (apply(Delta.mat, 2, function(delta.j) tapply(delta.j, factor(subject), FUN=sum)) != 0) + 0
}

#' Matrix with each element b_{kj} = I(\\sum_{i=1}^n \\delta_{ikj} \\ne 0)
#'
#' @export
fnGetBMat <- function(Delta.mat, condition) {
    (apply(Delta.mat, 2, function(delta.j) tapply(delta.j, factor(condition), FUN=sum)) != 0) + 0
}

#' Matrix with each element is \\Delta_{kj} = \\sum_{k^\\prime \\ne k} b_{kj}
#'
#' @export
fnGetDeltaSumMat <- function(B.mat) {
    #
    K <- nrow(B.mat)
    matrix(colSums(B.mat), nrow=1)[rep(1, K),] - B.mat
}

#' Calculate Mu.mat from the other parameters
#'
#' @export
fnGetMuMat <- function(r.vec, Alpha.mat, Theta.mat, Delta.mat, s.vec, i.vec, k.vec, X.mat=NA, Beta.mat=NA) {
    if(is.na(X.mat) && !is.na(Beta.mat)) stop("X.mat not NA but Beta.m provided")
    if(!is.na(X.mat) && is.na(Beta.mat)) stop("X.mat not NA but X.mat provided")

    J <- length(s.vec)
    K <- nrow(Theta.mat)
    n <- nrow(Alpha.mat)

    # Repeat the params as necessary to make NxJ matrices for Mu.mat
    R.mat.tmp <- matrix(r.vec, ncol=1)[, rep(1, J)]
    Alpha.mat.tmp <- Alpha.mat[i.vec,]
    Theta.mat.tmp <- Theta.mat[k.vec,]

    if(is.na(X.mat) && is.na(Beta.mat)) {
        Mu.mat <- exp(R.mat.tmp + Alpha.mat.tmp + Theta.mat.tmp)
    } else {
        Mu.mat <- exp(R.mat.tmp + Alpha.mat.tmp + Theta.mat.tmp + X.mat %*% Beta.mat)
    }
    Mu.mat[Delta.mat==0] <- 0

    return(Mu.mat)
}
