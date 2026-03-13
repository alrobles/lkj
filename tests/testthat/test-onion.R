# Tests for onion_cholesky, onion_corr, and ronion
library(testthat)
# ── onion_cholesky ────────────────────────────────────────────────────────────

test_that("onion_cholesky returns 1x1 identity for d=1", {
  L <- onion_cholesky(numeric(0), 1L, eta = 1)
  expect_true(is.matrix(L))
  expect_equal(dim(L), c(1L, 1L))
  expect_equal(L[1L, 1L], 1.0)
})

test_that("onion_cholesky returns valid 2x2 Cholesky factor", {
  L <- onion_cholesky(0.0, 2L, eta = 1)
  expect_true(is.matrix(L))
  expect_equal(dim(L), c(2L, 2L))
  expect_true(all(L[upper.tri(L)] == 0))
  expect_true(all(diag(L) > 0))
  R <- L %*% t(L)
  expect_equal(diag(R), c(1.0, 1.0), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

test_that("onion_cholesky returns lower-triangular matrix with positive diagonal", {
  d <- 4L
  m <- d * (d + 1L) / 2L - 2L   # = 8
  v <- rep(0.5, m)
  L <- onion_cholesky(v, d, eta = 1)
  expect_true(is.matrix(L))
  expect_equal(dim(L), c(d, d))
  expect_true(all(L[upper.tri(L)] == 0))
  expect_true(all(diag(L) > 0))
})

test_that("onion_cholesky produces a valid correlation matrix (R = L L')", {
  d <- 4L
  m <- d * (d + 1L) / 2L - 2L   # = 8
  v <- seq(-1, 1, length.out = m)
  L <- onion_cholesky(v, d, eta = 1)
  R <- L %*% t(L)
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
  expect_equal(R, t(R), tolerance = 1e-10)
  eigs <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("onion_cholesky handles d=3 correctly", {
  d <- 3L
  m <- d * (d + 1L) / 2L - 2L   # = 4
  v <- rep(0.0, m)
  L <- onion_cholesky(v, d, eta = 1)
  R <- L %*% t(L)
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
  eigs <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("onion_cholesky validates inputs", {
  expect_error(onion_cholesky(numeric(0), 0L),        "d must be >= 1")
  expect_error(onion_cholesky(numeric(0), 2L),        "v must have length")
  expect_error(onion_cholesky(c(0.1), 2L, eta = 0),  "eta must be > 0")
})

test_that("onion_cholesky with eta > 1 still produces valid correlation matrix", {
  d <- 5L
  m <- d * (d + 1L) / 2L - 2L   # = 13
  set.seed(42L)
  v <- rlogis(m)
  L <- onion_cholesky(v, d, eta = 4)
  R <- L %*% t(L)
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
  eigs <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("onion_cholesky is deterministic (same v gives same L)", {
  d <- 4L
  m <- d * (d + 1L) / 2L - 2L
  set.seed(7L)
  v <- rlogis(m)
  L1 <- onion_cholesky(v, d, eta = 2)
  L2 <- onion_cholesky(v, d, eta = 2)
  expect_equal(L1, L2)
})

# ── onion_corr ────────────────────────────────────────────────────────────────

test_that("onion_corr returns 1x1 identity for d=1", {
  R <- onion_corr(1L, eta = 1)
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(1L, 1L))
  expect_equal(R[1L, 1L], 1.0)
})

test_that("onion_corr returns a valid 2x2 correlation matrix", {
  set.seed(1L)
  R <- onion_corr(2L, eta = 1)
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(2L, 2L))
  expect_equal(diag(R), c(1.0, 1.0), tolerance = 1e-10)
  expect_equal(R, t(R), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

test_that("onion_corr returns a valid d x d correlation matrix", {
  set.seed(99L)
  d <- 3L
  R <- onion_corr(d, eta = 2)
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(d, d))
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
  expect_equal(R, t(R), tolerance = 1e-4)
  eigs <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("onion_corr with eta = 1 (uniform) still yields PD matrix", {
  set.seed(11L)
  R <- onion_corr(6L, eta = 1)
  eigs <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigs > 0))
  expect_equal(diag(R), rep(1.0, 6L), tolerance = 1e-10)
})

test_that("onion_corr off-diagonal entries are in (-1, 1)", {
  set.seed(55L)
  R <- onion_corr(4L, eta = 1)
  off_diag <- R[row(R) != col(R)]
  expect_true(all(off_diag > -1) && all(off_diag < 1))
})

test_that("onion_corr produces identical results with the same seed", {
  set.seed(21L)
  R1 <- onion_corr(4L, eta = 1)
  set.seed(21L)
  R2 <- onion_corr(4L, eta = 1)
  expect_equal(R1, R2)
})

# ── ronion ────────────────────────────────────────────────────────────────────

test_that("ronion returns a d x d matrix for n=1", {
  set.seed(123L)
  R <- ronion(1L, d = 3L, eta = 1)
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(3L, 3L))
  expect_equal(diag(R), rep(1.0, 3L), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

test_that("ronion returns a 3D array for n > 1", {
  set.seed(7L)
  arr <- ronion(5L, d = 4L, eta = 2)
  expect_equal(dim(arr), c(4L, 4L, 5L))
  for (i in seq_len(5L)) {
    R_i <- arr[, , i]
    expect_equal(diag(R_i), rep(1.0, 4L), tolerance = 1e-10)
    expect_true(all(eigen(R_i, only.values = TRUE)$values > 0))
  }
})

test_that("ronion with eta > 1 returns valid matrices", {
  set.seed(42L)
  arr <- ronion(3L, d = 3L, eta = 5)
  expect_equal(dim(arr), c(3L, 3L, 3L))
  for (i in seq_len(3L)) {
    R_i <- arr[, , i]
    expect_equal(diag(R_i), rep(1.0, 3L), tolerance = 1e-10)
    eigs <- eigen(R_i, only.values = TRUE)$values
    expect_true(all(eigs > 0))
  }
})
