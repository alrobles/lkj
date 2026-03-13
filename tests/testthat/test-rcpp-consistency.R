# Tests verifying that the C++ (Rcpp) implementations produce results
# consistent with the pure-R reference implementations.
library(testthat)

# ── Pure-R reference implementations (kept for verification) ──────────────────

.cvine_cholesky_r <- function(v, d, eta = 1) {
  stopifnot(d >= 1, eta > 0, is.numeric(v))
  if (d == 1L) return(matrix(1, 1, 1))
  m_needed <- d * (d - 1L) / 2L
  if (length(v) != m_needed)
    stop(sprintf("v must have length %d for d=%d (got %d).", m_needed, d, length(v)))

  to_unit <- function(x, eps = 1e-12) pmin(pmax(stats::plogis(x), eps), 1 - eps)
  p   <- matrix(0, d, d)
  idx <- 1L
  for (k in 1:(d - 1L)) {
    phi_k <- eta + (d - k - 1) / 2
    for (ell in (k + 1L):d) {
      u  <- to_unit(v[idx]); idx <- idx + 1L
      val <- 2 * stats::qbeta(u, phi_k, phi_k) - 1
      p[k, ell] <- max(min(val, 1 - 1e-15), -1 + 1e-15)
    }
  }
  L <- diag(d)
  for (j in 2:d) {
    r <- numeric(j - 1L)
    r[1] <- p[1, j]
    if (j >= 3L) {
      for (i in 2:(j - 1L)) {
        rij <- p[i, j]
        for (m in (i - 1L):1L) {
          rim <- p[m, i]; rjm <- p[m, j]
          rij <- rij * sqrt(pmax(0, (1 - rim^2) * (1 - rjm^2))) + rim * rjm
        }
        r[i] <- rij
      }
    }
    g <- forwardsolve(L[1:(j-1L), 1:(j-1L), drop=FALSE], r)
    L[j, 1:(j-1L)] <- g
    L[j, j] <- sqrt(pmax(0, 1 - sum(g^2)))
  }
  L
}

.onion_cholesky_r <- function(v, d, eta = 1) {
  stopifnot(d >= 1, eta > 0, is.numeric(v))
  d <- as.integer(d)
  if (d == 1L) return(matrix(1, 1, 1))
  m_needed <- if (d == 2L) 1L else d * (d + 1L) / 2L - 2L
  if (length(v) != m_needed)
    stop(sprintf("v must have length %d for d=%d (got %d).", m_needed, d, length(v)))

  to_unit <- function(x, eps = 1e-12) pmin(pmax(stats::plogis(x), eps), 1 - eps)
  L   <- matrix(0, d, d)
  L[1L, 1L] <- 1
  idx <- 1L
  beta <- eta + (d - 2) / 2

  u1  <- to_unit(v[idx]); idx <- idx + 1L
  r12 <- 2 * stats::qbeta(u1, beta, beta) - 1
  r12 <- max(min(r12, 1 - 1e-15), -1 + 1e-15)
  L[2L, 1L] <- r12
  L[2L, 2L] <- sqrt(max(0, 1 - r12^2))
  if (d == 2L) return(L)

  for (n in 2L:(d - 1L)) {  # n = dimension of existing submatrix
    beta  <- beta - 0.5
    k     <- n              # dimension of existing submatrix (same as n)
    row_j <- n + 1L
    u_y <- to_unit(v[idx]); idx <- idx + 1L
    y   <- max(min(stats::qbeta(u_y, k / 2, beta), 1 - 1e-15), 1e-15)
    z_raw <- stats::qnorm(to_unit(v[idx:(idx + k - 1L)]))
    idx   <- idx + k
    nrm   <- sqrt(sum(z_raw^2))
    u_vec <- if (nrm < 1e-12) c(1, rep(0, k - 1)) else z_raw / nrm
    L[row_j, 1L:k]  <- sqrt(y) * u_vec
    L[row_j, row_j] <- sqrt(max(0, 1 - y))
  }
  L
}

# ── cvine_cholesky: C++ vs R ───────────────────────────────────────────────────

test_that("cvine_cholesky C++ matches R for d=2", {
  v <- c(0.3)
  L_r   <- .cvine_cholesky_r(v, 2L, eta = 1)
  L_cpp <- cvine_cholesky_cpp(v, 2L, eta = 1)
  expect_equal(L_cpp, L_r, tolerance = 1e-10)
})

test_that("cvine_cholesky C++ matches R for d=4, eta=1", {
  set.seed(11L)
  d <- 4L
  v <- stats::rlogis(d * (d - 1L) / 2L)
  L_r   <- .cvine_cholesky_r(v, d, eta = 1)
  L_cpp <- cvine_cholesky_cpp(v, d, eta = 1)
  expect_equal(L_cpp, L_r, tolerance = 1e-10)
})

test_that("cvine_cholesky C++ matches R for d=5, eta=3", {
  set.seed(77L)
  d <- 5L
  v <- stats::rlogis(d * (d - 1L) / 2L)
  L_r   <- .cvine_cholesky_r(v, d, eta = 3)
  L_cpp <- cvine_cholesky_cpp(v, d, eta = 3)
  expect_equal(L_cpp, L_r, tolerance = 1e-9)
})

test_that("cvine_cholesky C++ matches R for d=1", {
  L_r   <- .cvine_cholesky_r(numeric(0), 1L)
  L_cpp <- cvine_cholesky_cpp(numeric(0), 1L)
  expect_equal(L_cpp, L_r, tolerance = 1e-12)
})

test_that("cvine_cholesky C++ produces valid correlation matrix for d=6", {
  set.seed(33L)
  d <- 6L
  v <- stats::rlogis(d * (d - 1L) / 2L)
  L <- cvine_cholesky_cpp(v, d, eta = 2)
  R <- L %*% t(L)
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

# ── onion_cholesky: C++ vs R ───────────────────────────────────────────────────

test_that("onion_cholesky C++ matches R for d=2", {
  v <- c(0.5)
  L_r   <- .onion_cholesky_r(v, 2L, eta = 1)
  L_cpp <- onion_cholesky_cpp(v, 2L, eta = 1)
  expect_equal(L_cpp, L_r, tolerance = 1e-10)
})

test_that("onion_cholesky C++ matches R for d=4, eta=1", {
  set.seed(11L)
  d <- 4L
  v <- stats::rlogis(d * (d + 1L) / 2L - 2L)
  L_r   <- .onion_cholesky_r(v, d, eta = 1)
  L_cpp <- onion_cholesky_cpp(v, d, eta = 1)
  expect_equal(L_cpp, L_r, tolerance = 1e-10)
})

test_that("onion_cholesky C++ matches R for d=5, eta=3", {
  set.seed(77L)
  d <- 5L
  v <- stats::rlogis(d * (d + 1L) / 2L - 2L)
  L_r   <- .onion_cholesky_r(v, d, eta = 3)
  L_cpp <- onion_cholesky_cpp(v, d, eta = 3)
  expect_equal(L_cpp, L_r, tolerance = 1e-9)
})

test_that("onion_cholesky C++ matches R for d=1", {
  L_r   <- .onion_cholesky_r(numeric(0), 1L)
  L_cpp <- onion_cholesky_cpp(numeric(0), 1L)
  expect_equal(L_cpp, L_r, tolerance = 1e-12)
})

test_that("onion_cholesky C++ produces valid correlation matrix for d=6", {
  set.seed(55L)
  d <- 6L
  v <- stats::rlogis(d * (d + 1L) / 2L - 2L)
  L <- onion_cholesky_cpp(v, d, eta = 2)
  R <- L %*% t(L)
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})
