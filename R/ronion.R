#' Generate multiple correlation matrices via the onion method
#'
#' @param n Number of matrices to generate.
#' @param d Dimension of the matrices.
#' @param eta Positive shape parameter.
#'
#' @return If \code{n = 1}, a \code{d x d} correlation matrix.
#'   If \code{n > 1}, an array of dimension \code{c(d, d, n)}.
#'
#' @examples
#' R_list <- ronion(5, d = 3, eta = 1)
#'
#' @export
ronion <- function(n, d, eta = 1) {
  stopifnot(
    n == as.integer(n), n >= 1,
    d == as.integer(d), d >= 1,
    is.numeric(eta), length(eta) == 1L, eta > 0
  )

  if (n == 1) {
    return(onion_corr(d, eta))
  }

  out <- array(0, dim = c(d, d, n))
  for (i in seq_len(n)) {
    out[, , i] <- onion_corr(d, eta)
  }
  out
}
