# Tests for the unified lkj() and rlkj() API
library(testthat)

# ── lkj() ─────────────────────────────────────────────────────────────────────

test_that("lkj() cvine returns a valid correlation matrix", {
  d <- 4L
  v <- c(-0.5, 0.3, 1.2, -0.8, 0.1, 0.6)
  R <- lkj(v, d = d, eta = 1, method = "cvine")
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(d, d))
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
  expect_equal(R, t(R), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

test_that("lkj() onion returns a valid correlation matrix", {
  d <- 4L
  v <- rep(0, d * (d + 1L) / 2L - 2L)
  R <- lkj(v, d = d, eta = 1, method = "onion")
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(d, d))
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
  expect_equal(R, t(R), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

test_that("lkj() output='L' returns valid lower-triangular Cholesky factor", {
  d <- 3L
  v <- rep(0, d * (d - 1L) / 2L)
  L <- lkj(v, d = d, eta = 1, method = "cvine", output = "L")
  expect_true(is.matrix(L))
  expect_equal(dim(L), c(d, d))
  expect_true(all(L[upper.tri(L)] == 0))
  expect_true(all(diag(L) > 0))
  R <- L %*% t(L)
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
})

test_that("lkj() output='R' matches L %*% t(L) for cvine", {
  d <- 3L
  v <- c(0.1, -0.3, 0.5)
  L <- lkj(v, d = d, eta = 1, method = "cvine", output = "L")
  R_from_lkj <- lkj(v, d = d, eta = 1, method = "cvine", output = "R")
  expect_equal(R_from_lkj, L %*% t(L), tolerance = 1e-12)
})

test_that("lkj() output='R' matches L %*% t(L) for onion", {
  d <- 3L
  v <- rep(0.2, d * (d + 1L) / 2L - 2L)
  L <- lkj(v, d = d, eta = 1, method = "onion", output = "L")
  R_from_lkj <- lkj(v, d = d, eta = 1, method = "onion", output = "R")
  expect_equal(R_from_lkj, L %*% t(L), tolerance = 1e-12)
})

test_that("lkj() handles d=1 for both methods", {
  R_cv <- lkj(numeric(0), d = 1L, eta = 1, method = "cvine")
  R_on <- lkj(numeric(0), d = 1L, eta = 1, method = "onion")
  expect_equal(R_cv, matrix(1, 1, 1))
  expect_equal(R_on, matrix(1, 1, 1))
})

test_that("lkj() with eta > 1 still produces valid matrices", {
  d <- 5L
  v_cv <- rep(0, d * (d - 1L) / 2L)
  v_on <- rep(0, d * (d + 1L) / 2L - 2L)
  R_cv <- lkj(v_cv, d = d, eta = 5, method = "cvine")
  R_on <- lkj(v_on, d = d, eta = 5, method = "onion")
  expect_equal(diag(R_cv), rep(1.0, d), tolerance = 1e-10)
  expect_equal(diag(R_on), rep(1.0, d), tolerance = 1e-10)
  expect_true(all(eigen(R_cv, only.values = TRUE)$values > 0))
  expect_true(all(eigen(R_on, only.values = TRUE)$values > 0))
})

# ── rlkj() ────────────────────────────────────────────────────────────────────

test_that("rlkj() cvine n=1 returns a valid correlation matrix", {
  set.seed(42L)
  R <- rlkj(1L, d = 3L, eta = 1, method = "cvine")
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(3L, 3L))
  expect_equal(diag(R), rep(1.0, 3L), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

test_that("rlkj() onion n=1 returns a valid correlation matrix", {
  set.seed(42L)
  R <- rlkj(1L, d = 3L, eta = 1, method = "onion")
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(3L, 3L))
  expect_equal(diag(R), rep(1.0, 3L), tolerance = 1e-10)
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

test_that("rlkj() cvine n>1 returns a 3D array of valid correlation matrices", {
  set.seed(7L)
  arr <- rlkj(5L, d = 4L, eta = 2, method = "cvine")
  expect_equal(dim(arr), c(4L, 4L, 5L))
  for (i in seq_len(5L)) {
    R_i <- arr[, , i]
    expect_equal(diag(R_i), rep(1.0, 4L), tolerance = 1e-10)
    expect_true(all(eigen(R_i, only.values = TRUE)$values > 0))
  }
})

test_that("rlkj() onion n>1 returns a 3D array of valid correlation matrices", {
  set.seed(7L)
  arr <- rlkj(5L, d = 4L, eta = 2, method = "onion")
  expect_equal(dim(arr), c(4L, 4L, 5L))
  for (i in seq_len(5L)) {
    R_i <- arr[, , i]
    expect_equal(diag(R_i), rep(1.0, 4L), tolerance = 1e-10)
    expect_true(all(eigen(R_i, only.values = TRUE)$values > 0))
  }
})

test_that("rlkj() output='L' returns valid Cholesky factors", {
  set.seed(99L)
  arr <- rlkj(3L, d = 3L, eta = 1, method = "cvine", output = "L")
  expect_equal(dim(arr), c(3L, 3L, 3L))
  for (i in seq_len(3L)) {
    L_i <- arr[, , i]
    expect_true(all(L_i[upper.tri(L_i)] == 0))
    expect_true(all(diag(L_i) > 0))
    R_i <- L_i %*% t(L_i)
    expect_equal(diag(R_i), rep(1.0, 3L), tolerance = 1e-10)
  }
})

test_that("rlkj() handles d=1", {
  R <- rlkj(1L, d = 1L, eta = 1, method = "cvine")
  expect_equal(R, matrix(1, 1, 1))
  R2 <- rlkj(1L, d = 1L, eta = 1, method = "onion")
  expect_equal(R2, matrix(1, 1, 1))
})

test_that("rlkj() default method is cvine", {
  set.seed(123L)
  R1 <- rlkj(1L, d = 3L, eta = 1)
  set.seed(123L)
  R2 <- rlkj(1L, d = 3L, eta = 1, method = "cvine")
  expect_equal(R1, R2)
})
