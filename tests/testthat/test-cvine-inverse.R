library(testthat)

# ── cvine_inverse: basic properties ───────────────────────────────────────────

test_that("cvine_inverse returns numeric(0) for d=1", {
  v <- cvine_inverse(matrix(1, 1, 1), eta = 1)
  expect_true(is.numeric(v))
  expect_equal(length(v), 0L)
})

test_that("cvine_inverse returns a vector of length d*(d-1)/2", {
  for (d in c(2L, 3L, 4L, 5L, 6L)) {
    set.seed(d)
    v_in <- stats::rlogis(d * (d - 1L) / 2L)
    L    <- cvine_cholesky(v_in, d, eta = 1)
    R    <- L %*% t(L)
    v_out <- cvine_inverse(R, eta = 1)
    expect_equal(length(v_out), d * (d - 1L) / 2L,
                 label = paste0("d=", d))
  }
})

# ── cvine_inverse: round-trip v -> R -> v ─────────────────────────────────────

test_that("cvine_inverse recovers the original v for d=2, eta=1", {
  v <- c(0.5)
  L <- cvine_cholesky(v, 2L, eta = 1)
  R <- L %*% t(L)
  v_back <- cvine_inverse(R, eta = 1)
  expect_equal(v_back, v, tolerance = 1e-10)
})

test_that("cvine_inverse recovers the original v for d=4, eta=1", {
  v <- c(-0.5, 0.3, 1.2, -0.8, 0.1, 0.6)
  L <- cvine_cholesky(v, 4L, eta = 1)
  R <- L %*% t(L)
  v_back <- cvine_inverse(R, eta = 1)
  expect_equal(v_back, v, tolerance = 1e-10)
})

test_that("cvine_inverse recovers the original v for d=5, eta=2", {
  set.seed(42L)
  d <- 5L
  v <- stats::rlogis(d * (d - 1L) / 2L)
  L <- cvine_cholesky(v, d, eta = 2)
  R <- L %*% t(L)
  v_back <- cvine_inverse(R, eta = 2)
  expect_equal(v_back, v, tolerance = 1e-9)
})

test_that("cvine_inverse recovers the original v for d=6, eta=3", {
  set.seed(77L)
  d <- 6L
  v <- stats::rlogis(d * (d - 1L) / 2L)
  L <- cvine_cholesky(v, d, eta = 3)
  R <- L %*% t(L)
  v_back <- cvine_inverse(R, eta = 3)
  expect_equal(v_back, v, tolerance = 1e-9)
})

test_that("cvine_inverse handles all-zero v (identity-like matrix) for d=3", {
  d <- 3L
  v <- rep(0.0, d * (d - 1L) / 2L)
  L <- cvine_cholesky(v, d, eta = 1)
  R <- L %*% t(L)
  v_back <- cvine_inverse(R, eta = 1)
  expect_equal(v_back, v, tolerance = 1e-10)
})

# ── cvine_inverse: round-trip R -> v -> R ─────────────────────────────────────

test_that("cvine_inverse round-trip R->v->R for a randomly generated matrix", {
  set.seed(123L)
  R_in <- rlkj_cvine(1L, d = 4L, eta = 1)
  v    <- cvine_inverse(R_in, eta = 1)
  L    <- cvine_cholesky(v, 4L, eta = 1)
  R_out <- L %*% t(L)
  expect_equal(R_out, R_in, tolerance = 1e-10)
})

test_that("cvine_inverse round-trip R->v->R for d=5 with eta=2", {
  set.seed(55L)
  d    <- 5L
  R_in <- rlkj_cvine(1L, d = d, eta = 2)
  v    <- cvine_inverse(R_in, eta = 2)
  L    <- cvine_cholesky(v, d, eta = 2)
  R_out <- L %*% t(L)
  expect_equal(R_out, R_in, tolerance = 1e-10)
})

# ── cvine_inverse: input validation ───────────────────────────────────────────

test_that("cvine_inverse rejects non-matrix input", {
  expect_error(cvine_inverse(c(1, 0.5, 0.5, 1)), "R must be a matrix")
})

test_that("cvine_inverse rejects non-square matrix", {
  expect_error(cvine_inverse(matrix(1:6, 2, 3)), "R must be square")
})

test_that("cvine_inverse rejects non-unit diagonal", {
  R_bad <- matrix(c(2, 0.5, 0.5, 2), 2, 2)
  expect_error(cvine_inverse(R_bad, eta = 1), "unit diagonal")
})

test_that("cvine_inverse rejects non-symmetric matrix", {
  R_bad <- matrix(c(1, 0.5, 0.3, 1), 2, 2)
  expect_error(cvine_inverse(R_bad, eta = 1), "symmetric")
})

test_that("cvine_inverse rejects non-positive-definite matrix", {
  # A near-unit correlation (valid PD) should work without error
  R_ok <- matrix(c(1, 0.999, 0.999, 1), 2, 2)
  expect_no_error(cvine_inverse(R_ok, eta = 1))
  # A singular matrix (eigenvalue = 0) is not positive definite
  R_npd <- matrix(c(1, -1, -1, 1), 2, 2)
  expect_error(cvine_inverse(R_npd, eta = 1), "positive definite")
})

test_that("cvine_inverse rejects eta <= 0", {
  R <- diag(3)
  R[1, 2] <- R[2, 1] <- 0.3
  R[1, 3] <- R[3, 1] <- 0.2
  R[2, 3] <- R[3, 2] <- 0.1
  expect_error(cvine_inverse(R, eta = 0),  "eta must be > 0")
  expect_error(cvine_inverse(R, eta = -1), "eta must be > 0")
})

# ── cvine_inverse: C++ implementation consistency ────────────────────────────

test_that("cvine_inverse_cpp and cvine_inverse agree", {
  set.seed(99L)
  d <- 4L
  v <- stats::rlogis(d * (d - 1L) / 2L)
  L <- cvine_cholesky(v, d, eta = 1)
  R <- L %*% t(L)
  v_r   <- cvine_inverse(R, eta = 1)
  v_cpp <- cvine_inverse_cpp(R, d, eta = 1.0)
  expect_equal(v_r, v_cpp, tolerance = 1e-14)
})

# ── cvine_inverse: determinism ────────────────────────────────────────────────

test_that("cvine_inverse is deterministic", {
  set.seed(7L)
  d <- 4L
  v <- stats::rlogis(d * (d - 1L) / 2L)
  L <- cvine_cholesky(v, d, eta = 1)
  R <- L %*% t(L)
  v1 <- cvine_inverse(R, eta = 1)
  v2 <- cvine_inverse(R, eta = 1)
  expect_equal(v1, v2)
})
