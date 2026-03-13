# Build Cholesky factor L directly from C-vine partial correlations (LKJ).
# - Uses eta > 0 to set per-level Beta shapes: phi_k = eta + (d - k - 1)/2
# - v supplies d*(d-1)/2 reals, one per C-vine edge (in level-major order)
# - Returns lower-triangular L with diag(L) > 0 such that R = L %*% t(L)
#
# Reference: Lewandowski, Kurowicka & Joe (2009), Sec. 2.4 (C-vine) and Eq. (2).

cvine_cholesky <- function(v, d, eta = 1) {
  stopifnot(d >= 1, eta > 0, is.numeric(v))

  if (d == 1L) return(matrix(1, 1, 1))

  # Need exactly one real per edge: m = d*(d-1)/2
  m_needed <- d * (d - 1L) / 2L
  if (length(v) != m_needed) {
    stop(sprintf("v must have length %d for d=%d (got %d).",
                 m_needed, d, length(v)))
  }

  # Stable map to (0,1)
  to_unit <- function(x, eps = 1e-12) pmin(pmax(plogis(x), eps), 1 - eps)

  # Allocate partial correlation matrix p, where
  # p[k, l] = rho_{k,l | 1:(k-1)} for 1 <= k < l <= d in the C-vine order 1..d
  p <- matrix(0, d, d)

  # Fill p level by level with symmetric Beta on (-1,1):
  # For level k, shape is phi_k = eta + (d - k - 1)/2.
  idx <- 1L
  for (k in 1:(d - 1L)) {
    phi_k <- eta + (d - k - 1) / 2
    for (ell in (k + 1L):d) {
      u  <- to_unit(v[idx]); idx <- idx + 1L
      val <- 2 * qbeta(u, phi_k, phi_k) - 1
      # clamp slightly to avoid exactly +/-1
      p[k, ell] <- max(min(val, 1 - 1e-15), -1 + 1e-15)
    }
  }

  # Initialize L
  L <- diag(d)

  # j = 2..d: build each new row of L
  for (j in 2:d) {
    # Compute the correlation vector r with the previous variables
    r <- numeric(j - 1L)

    # r[1] = rho_{1j} = p[1, j]
    r[1] <- p[1, j]

    if (j >= 3L) {
      for (i in 2:(j - 1L)) {
        # Start from highest-order partial rho_{i,j | 1:(i-1)} = p[i, j]
        rij <- p[i, j]
        # Remove conditioning variables m = i-1 down to 1 via Yule–Kendall recursion
        # rho_{ij | L \ {m}} = rho_{ij|L} * sqrt((1-r_im^2)(1-r_jm^2)) + r_im * r_jm,
        # where r_im = p[m, i], r_jm = p[m, j] (both are partials with 1:(m-1) conditioned).
        for (m in (i - 1L):1L) {
          rim <- p[m, i]
          rjm <- p[m, j]
          rij <- rij * sqrt(pmax(0, (1 - rim^2) * (1 - rjm^2))) + rim * rjm
        }
        r[i] <- rij
      }
    }

    # Solve for the next row of L: L[1:(j-1), 1:(j-1)] %*% g = r
    # Use forwardsolve because the block is lower triangular with positive diagonal
    g <- forwardsolve(L[1:(j - 1L), 1:(j - 1L), drop = FALSE], r, upper.tri = FALSE, transpose = FALSE)

    # Append row j
    L[j, 1:(j - 1L)] <- g
    L[j, j] <- sqrt(pmax(0, 1 - sum(g^2)))
  }

  L
}
