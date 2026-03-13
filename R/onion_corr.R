#' Generate a correlation matrix via the extended onion method
#'
#' Implements the onion method of Ghosh & Henderson (2003) as extended by
#' Lewandowski, Kurowicka & Joe (2009). The resulting correlation matrix
#' has density proportional to \eqn{\det(\mathbf{R})^{\eta-1}}.
#'
#' @param d Dimension of the matrix (>= 1).
#' @param eta Positive shape parameter (default 1 gives uniform distribution).
#'
#' @return A \code{d x d} positive definite correlation matrix.
#'
#' @references
#' Lewandowski, D., Kurowicka, D., & Joe, H. (2009).
#' Generating random correlation matrices based on vines and extended onion method.
#' Journal of Multivariate Analysis, 100(9), 1989–2001.
#'
#' @examples
#' R <- onion_corr(d = 4, eta = 2)
#' print(R)
#'
#' @export
onion_corr <- function(d, eta = 1) {
  stopifnot(
    d == as.integer(d), d >= 1L,
    is.numeric(eta), length(eta) == 1L, eta > 0
  )

  if (d == 1L) {
    return(matrix(1, 1, 1))
  }

  # Step 1: initialise with 2x2 matrix
  beta <- eta + (d - 2) / 2
  r12 <- 2 * stats::rbeta(1, beta, beta) - 1
  R <- matrix(c(1, r12, r12, 1), 2, 2)

  if (d == 2L) {
    return(R)
  }

  # Step 2: sequentially add rows/columns
  for (n in 2:(d - 1L)) {
    beta <- beta - 0.5                 # a) decrease beta
    k <- n                             # current dimension

    # b) generate y ~ Beta(k/2, beta)
    y <- stats::rbeta(1, k / 2, beta)

    # c) generate uniform vector on k‑sphere
    u <- stats::rnorm(k)
    u <- u / sqrt(sum(u^2))

    # d) w = sqrt(y) * u
    w <- sqrt(y) * u

    # e) Cholesky factor A such that A A^T = R (lower triangular)
    A <- t(chol(R))                    # lower triangular

    # f) z = A %*% w
    z <- as.vector(A %*% w)

    # g) expand R
    R <- rbind(cbind(R, z), c(z, 1))
  }

  R
}
