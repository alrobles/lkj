library(testthat)
test_that("cvine_cholesky returns a valid Cholesky factor for d=1", {
  L <- cvine_cholesky(numeric(0), 1L, eta = 1)
  expect_true(is.matrix(L))
  expect_equal(dim(L), c(1L, 1L))
  expect_equal(L[1, 1], 1.0)
})

test_that("cvine_cholesky returns a lower-triangular matrix with positive diagonal", {
  d <- 4L
  v <- rep(0.5, d * (d - 1L) / 2L)
  L <- cvine_cholesky(v, d, eta = 1)
  expect_true(is.matrix(L))
  expect_equal(dim(L), c(d, d))
  # Lower triangular: upper triangle is zero
  expect_true(all(L[upper.tri(L)] == 0))
  # Positive diagonal
  expect_true(all(diag(L) > 0))
})

test_that("cvine_cholesky produces a valid correlation matrix (R = L L')", {
  d <- 4L
  v <- c(-0.5, 0.3, 1.2, -0.8, 0.1, 0.6)
  L <- cvine_cholesky(v, d, eta = 1)
  R <- L %*% t(L)
  # Diagonal of R must be 1
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
  # R must be symmetric
  expect_equal(R, t(R), tolerance = 1e-10)
  # R must be positive definite
  eigs <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("cvine_cholesky handles d=2 correctly", {
  d <- 2L
  v <- 0.0  # logistic(0) = 0.5, qbeta(0.5, ...) = 0.5 → val = 0
  L <- cvine_cholesky(v, d, eta = 1)
  R <- L %*% t(L)
  expect_equal(diag(R), c(1.0, 1.0), tolerance = 1e-10)
  eigs <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

# test_that("cvine_cholesky validates inputs", {
#   expect_error(cvine_cholesky(numeric(0), 0L),  "d must be >= 1")
#   expect_error(cvine_cholesky(numeric(0), 2L),  "v must have length")
#   expect_error(cvine_cholesky(c(0.1, 0.2), 2L, eta = 0),  "eta must be > 0")
# })

test_that("cvine_cholesky eta > 1 still produces valid correlation matrix", {
  d <- 5L
  m <- d * (d - 1L) / 2L
  set.seed(42L)
  v <- rlogis(m)
  L <- cvine_cholesky(v, d, eta = 5)
  R <- L %*% t(L)
  expect_equal(diag(R), rep(1.0, d), tolerance = 1e-10)
  eigs <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("rlkj_cvine returns a valid d x d correlation matrix for n=1", {
  set.seed(123L)
  R <- rlkj_cvine(1L, d = 3L, eta = 2)
  expect_true(is.matrix(R))
  expect_equal(dim(R), c(3L, 3L))
  expect_equal(diag(R), rep(1.0, 3L), tolerance = 1e-10)
  eigs <- eigen(R, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("rlkj_cvine returns a 3D array for n > 1", {
  set.seed(7L)
  arr <- rlkj_cvine(5L, d = 3L, eta = 1)
  expect_equal(dim(arr), c(3L, 3L, 5L))
  for (i in seq_len(5L)) {
    R_i <- arr[, , i]
    expect_equal(diag(R_i), rep(1.0, 3L), tolerance = 1e-10)
    expect_true(all(eigen(R_i, only.values = TRUE)$values > 0))
  }
})

