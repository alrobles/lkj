#' Generate random correlation matrices from the LKJ distribution via C-vine
#'
#' Draws \code{n} random correlation matrices of dimension \code{d} whose
#' density is proportional to \eqn{\det(R)^{\eta-1}}. The method uses the
#' C-vine construction implemented in C++ via \code{\link{cvine_cholesky}},
#' with independent logistic draws for the underlying unconstrained reals.
#'
#' @param n Number of matrices to generate.
#' @param d Dimension of the correlation matrices.
#' @param eta Positive shape parameter (default 1 gives the uniform
#'   distribution over correlation matrices).
#'
#' @return If \code{n = 1}, a \code{d x d} correlation matrix.
#'   If \code{n > 1}, an array of dimension \code{c(d, d, n)}.
#'
#' @examples
#' R <- rlkj_cvine(1, d = 4, eta = 2)
#' R_arr <- rlkj_cvine(10, d = 3, eta = 1)
#'
#' @export
rlkj_cvine <- function(n, d, eta = 1) {
  stopifnot(
    n == as.integer(n), n >= 1,
    d == as.integer(d), d >= 1,
    is.numeric(eta), length(eta) == 1L, eta > 0
  )

  m <- d * (d - 1L) / 2L
  if (n == 1L) {
    v <- stats::rlogis(m)
    L <- cvine_cholesky(v, d, eta)
    return(L %*% t(L))
  }

  # For n > 1, return a 3D array of correlation matrices
  out <- array(0, dim = c(d, d, n))
  for (i in seq_len(n)) {
    v <- stats::rlogis(m)
    L <- cvine_cholesky(v, d, eta)
    out[, , i] <- L %*% t(L)
  }
  out
}
