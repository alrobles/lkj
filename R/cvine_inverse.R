#' Inverse of the C-vine LKJ mapping: correlation matrix to unconstrained reals
#'
#' Maps a correlation matrix \eqn{R} back to the vector of unconstrained real
#' parameters \eqn{v} such that \code{cvine_cholesky(v, d, eta)} produces the
#' Cholesky factor \eqn{L} satisfying \eqn{R = L L^\top}.
#'
#' The inverse recovers the table of partial correlations from the off-diagonal
#' entries of \eqn{R} using the inverse C-vine recursion, then applies the
#' inverse Beta CDF and logit to map each partial correlation to an
#' unconstrained real.
#'
#' @param R A \code{d x d} numeric matrix that is a valid correlation matrix:
#'   symmetric, positive-definite, and with unit diagonal.
#' @param eta Positive numeric. LKJ shape parameter. Must match the value of
#'   \code{eta} used in the forward mapping to obtain an exact round-trip.
#'
#' @return A numeric vector of length \eqn{d(d-1)/2} of unconstrained reals.
#'   For \code{d = 1} an empty vector \code{numeric(0)} is returned.
#'
#' @references
#' Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random
#' correlation matrices based on vines and extended onion method.
#' \emph{Journal of Multivariate Analysis}, 100(9), 1989–2001.
#'
#' @examples
#' # Round-trip: v -> R -> v
#' v <- c(-0.5, 0.3, 1.2, -0.8, 0.1, 0.6)
#' L <- cvine_cholesky(v, d = 4L, eta = 1)
#' R <- L %*% t(L)
#' v_recovered <- cvine_inverse(R, eta = 1)
#' stopifnot(max(abs(v - v_recovered)) < 1e-10)
#'
#' @export
cvine_inverse <- function(R, eta = 1) {
  stopifnot(
    "R must be a matrix"        = is.matrix(R),
    "R must be square"          = nrow(R) == ncol(R),
    "R must be numeric"         = is.numeric(R),
    "eta must be > 0"           = is.numeric(eta) && length(eta) == 1L && eta > 0
  )

  d <- nrow(R)

  if (!isTRUE(all.equal(diag(R), rep(1.0, d), tolerance = 1e-8))) {
    stop("R must have unit diagonal (it must be a correlation matrix).")
  }
  if (!isTRUE(all.equal(R, t(R), tolerance = 1e-8))) {
    stop("R must be symmetric.")
  }
  eigs <- eigen(R, symmetric = TRUE, only.values = TRUE)$values
  if (any(eigs <= 0)) {
    stop("R must be positive definite (all eigenvalues must be > 0).")
  }

  if (d == 1L) return(numeric(0L))

  cvine_inverse_cpp(R, d, eta)
}
