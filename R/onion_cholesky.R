#' Deterministic onion-method Cholesky factor from real parameters
#'
#' Implements the deterministic mapping from a vector of real numbers to a
#' lower-triangular Cholesky factor \code{L} such that \code{R = L \%*\% t(L)}
#' is a valid correlation matrix with density proportional to
#' \eqn{\det(R)^{\eta-1}}.
#'
#' The mapping follows the sequential onion construction of Lewandowski,
#' Kurowicka & Joe (2009), replacing random draws with quantile functions
#' applied to logistic-transformed inputs (see \code{METHODS.md}, Section 5).
#' Parameters are consumed in order: the first real encodes \eqn{r_{12}};
#' at each subsequent step \eqn{k} one real encodes the radial component
#' \eqn{y \sim \mathrm{Beta}(k/2, \beta)} and \eqn{k} reals encode a
#' direction on the \eqn{k}-sphere via independent standard-normal projections.
#'
#' The required input length is \eqn{d(d+1)/2 - 2} for \eqn{d \ge 2} and
#' \eqn{0} for \eqn{d = 1}.
#'
#' @param v Numeric vector of reals. Length must equal \eqn{d(d+1)/2 - 2}
#'   for \eqn{d \ge 2} (0 for \eqn{d = 1}).
#' @param d Positive integer dimension of the target correlation matrix.
#' @param eta Positive shape parameter (default 1 gives the uniform
#'   distribution over correlation matrices).
#'
#' @return A \eqn{d \times d} lower-triangular matrix \code{L} with positive
#'   diagonal such that \code{R = L \%*\% t(L)} is a positive-definite
#'   correlation matrix.
#'
#' @references
#' Lewandowski, D., Kurowicka, D., & Joe, H. (2009).
#' Generating random correlation matrices based on vines and extended onion
#' method. \emph{Journal of Multivariate Analysis}, 100(9), 1989--2001.
#'
#' @seealso \code{\link{cvine_cholesky}} for the C-vine deterministic mapping,
#'   \code{\link{onion_corr}} for random generation via the onion method.
#'
#' @examples
#' ## d = 4: need d*(d+1)/2 - 2 = 8 parameters
#' v <- rep(0, 8)
#' L <- onion_cholesky(v, d = 4, eta = 1)
#' R <- L %*% t(L)
#' stopifnot(all(abs(diag(R) - 1) < 1e-10))
#'
#' @export
onion_cholesky <- function(v, d, eta = 1) {
  stopifnot(
    "d must be >= 1" = d >= 1,
    "eta must be > 0" = eta > 0,
    "v must be numeric" = is.numeric(v)
  )
  d <- as.integer(d)

  if (d == 1L) return(matrix(1, 1, 1))

  # Number of parameters: 1 (for r12) + sum_{k=2}^{d-1} (1 + k) = d*(d+1)/2 - 2
  m_needed <- if (d == 2L) 1L else d * (d + 1L) / 2L - 2L
  if (length(v) != m_needed) {
    stop(sprintf(
      "v must have length %d for d=%d (got %d).",
      m_needed, d, length(v)
    ))
  }

  # Stable map from R to (0, 1)
  to_unit <- function(x, eps = 1e-12) pmin(pmax(stats::plogis(x), eps), 1 - eps)

  # Initialise the lower-triangular Cholesky factor
  L   <- matrix(0, d, d)
  L[1L, 1L] <- 1
  idx <- 1L
  beta <- eta + (d - 2) / 2  # starting Beta shape

  # ── Step 1: r12 via symmetric Beta on (-1, 1) ──────────────────────────────
  u1  <- to_unit(v[idx]); idx <- idx + 1L
  r12 <- 2 * stats::qbeta(u1, beta, beta) - 1
  r12 <- max(min(r12, 1 - 1e-15), -1 + 1e-15)
  L[2L, 1L] <- r12
  L[2L, 2L] <- sqrt(max(0, 1 - r12^2))

  if (d == 2L) return(L)

  # ── Sequential steps: add dimensions 3, ..., d ─────────────────────────────
  # At step n (n = 2, ..., d-1) the current Cholesky block is n x n and we
  # append row n+1.  The mathematics (METHODS.md §3–5) shows that the new row
  # is w = sqrt(y) * u and the diagonal entry is sqrt(1 - y), where
  #   y  ~ Beta(k/2, beta)  (k = n)
  #   u  is the normalised direction (uniform on k-sphere)
  # so that z = A w is automatically in the correct ellipsoid and the new
  # Cholesky row equals w (since A^{-1} z = w).
  for (n in 2L:(d - 1L)) {
    beta  <- beta - 0.5
    k     <- n            # current dimension
    row_j <- n + 1L

    # ── Radial component y ────────────────────────────────────────────────────
    u_y <- to_unit(v[idx]); idx <- idx + 1L
    y   <- max(min(stats::qbeta(u_y, k / 2, beta), 1 - 1e-15), 1e-15)

    # ── Direction: k reals → k standard normals → unit vector ────────────────
    z_raw <- stats::qnorm(to_unit(v[idx:(idx + k - 1L)]))
    idx   <- idx + k
    nrm   <- sqrt(sum(z_raw^2))
    if (nrm < 1e-12) nrm <- 1  # guard against degenerate input
    u_vec <- z_raw / nrm

    # ── New row of L ──────────────────────────────────────────────────────────
    L[row_j, 1L:k]  <- sqrt(y) * u_vec
    L[row_j, row_j] <- sqrt(max(0, 1 - y))
  }

  L
}
