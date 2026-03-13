onion_cholesky <- function(v, d, eta = 1) {
  stopifnot(
    "d must be >= 1" = d >= 1,
    "eta must be > 0" = eta > 0,
    "v must be numeric" = is.numeric(v)
  )
  d <- as.integer(d)

  if (d == 1L) return(matrix(1, 1, 1))

  m_needed <- if (d == 2L) 1L else d * (d + 1L) / 2L - 2L
  if (length(v) != m_needed) {
    stop(sprintf(
      "v must have length %d for d=%d (got %d).",
      m_needed, d, length(v)
    ))
  }

  to_unit <- function(x, eps = 1e-12) pmin(pmax(stats::plogis(x), eps), 1 - eps)

  L   <- matrix(0, d, d)
  L[1L, 1L] <- 1
  idx <- 1L
  beta <- eta + (d - 2) / 2

  # Step 1: r12
  u1  <- to_unit(v[idx]); idx <- idx + 1L
  r12 <- 2 * stats::qbeta(u1, beta, beta) - 1
  r12 <- max(min(r12, 1 - 1e-15), -1 + 1e-15)
  L[2L, 1L] <- r12
  L[2L, 2L] <- sqrt(max(0, 1 - r12^2))

  if (d == 2L) return(L)

  # Sequential steps
  for (n in 2L:(d - 1L)) {
    beta  <- beta - 0.5
    k     <- n
    row_j <- n + 1L

    # Radial component y
    u_y <- to_unit(v[idx]); idx <- idx + 1L
    y   <- max(min(stats::qbeta(u_y, k / 2, beta), 1 - 1e-15), 1e-15)

    # Direction on the sphere
    z_raw <- stats::qnorm(to_unit(v[idx:(idx + k - 1L)]))
    idx   <- idx + k
    nrm   <- sqrt(sum(z_raw^2))

    if (nrm < 1e-12) {
      # Degenerate input – pick a fixed unit vector (first coordinate = 1)
      u_vec <- c(1, rep(0, k - 1))
    } else {
      u_vec <- z_raw / nrm
    }

    # New row of L
    L[row_j, 1L:k]  <- sqrt(y) * u_vec
    L[row_j, row_j] <- sqrt(max(0, 1 - y))
  }

  L
}
