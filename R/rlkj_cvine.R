#' Generate random correlation matrices from the LKJ distribution via C-vine
#'
#' Draws \code{n} random correlation matrices of dimension \code{d} whose density
#' is proportional to \eqn{\det(R)^{\eta-1}}. The method uses the C-vine
#' construction with independent logistic draws for the underlying partial
#' correlations.
#'
#' @param n Number of matrices to generate.
#' @param d Dimension of the matrices.
#' @param eta Positive shape parameter (default 1 gives uniform distribution).
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
  if (n == 1) {
    v <- rlogis(m)
    return(cvine_cholesky(v, d, eta) %*% t(cvine_cholesky(v, d, eta)))
  }

  # For n > 1, return a 3D array
  out <- array(0, dim = c(d, d, n))
  for (i in seq_len(n)) {
    v <- rlogis(m)
    L <- cvine_cholesky(v, d, eta)
    out[, , i] <- L %*% t(L)
  }
  out
}
