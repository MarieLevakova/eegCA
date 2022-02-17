##' Fitting reduced-rank regression with a specific rank
##'
##' Given a response matrix and a covariate matrix, this function fits reduced
##' rank regression for a specified rank. It reduces to singular value
##' decomposition if the covariate matrix is the identity matrix.
##'
##' @usage
##' rrr.fit(Y, X, nrank = 1, weight = NULL, coefSVD = FALSE)
##'
##' @param Y a matrix of response (n by q)
##' @param X a matrix of covariate (n by p)
##' @param nrank an integer specifying the desired rank
##' @param weight a square matrix of weight (q by q); The default is the
##'     identity matrix
##' @param coefSVD logical indicating the need for SVD for the coefficient matrix
##'     in the output; used in ssvd estimation
##' @return S3 \code{rrr} object, a list consisting of \item{coef}{coefficient
##'     of rrr} \item{coef.ls}{coefficient of least square} \item{fitted}{fitted
##'     value of rrr} \item{fitted.ls}{fitted value of least square}
##'     \item{A}{right singular matrix} \item{Ad}{a vector of sigular values}
##'     \item{rank}{rank of the fitted rrr}
##' @examples
##' Y <- matrix(rnorm(400), 100, 4)
##' X <- matrix(rnorm(800), 100, 8)
##' rfit <- rrr.fit(Y, X, nrank = 2)
##' coef(rfit)
##' @importFrom MASS ginv
##' @export
##'

rrr.fit <- function(Y, X, nrank = 1, coefSVD = FALSE){

  q <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  stopifnot(n == nrow(X))

  S_yx <- crossprod(Y, X)
  S_xx <- crossprod(X)

  ## FIXME: 0.01 is too arbitrary
  S_xx_inv <- tryCatch(
    ginv(S_xx),
    error = function(e)
      solve(S_xx + 0.01 * diag(p))
  )

  C_ls <- tcrossprod(S_xx_inv, S_yx)

  XC <- X %*% C_ls
  svdXC <- svd(XC, nrank, nrank)
  A <- svdXC$v[, 1:nrank]
  Ad <- (svdXC$d[1:nrank]) ^ 2
  AA <- tcrossprod(A)
  C_rr <- C_ls %*% AA

  ret <- list(
    coef = C_rr,
    coef.ls = C_ls,
    fitted = X %*% C_rr,
    fitted.ls = XC,
    A = A,
    Ad = Ad,
    rank = nrank
  )

  if (coefSVD) {
    coefSVD <- svd(C_rr, nrank, nrank)
    coefSVD$d <- coefSVD$d[1:nrank]
    coefSVD$u <- coefSVD$u[, 1:nrank, drop = FALSE]
    coefSVD$v <- coefSVD$v[, 1:nrank, drop = FALSE]
    ret <- c(ret, list(coefSVD = coefSVD))
  }

  ret
}
