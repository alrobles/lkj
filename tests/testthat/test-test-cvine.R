test_that("cvine_cholesky returns a valid Cholesky factor", {
  d <- 4
  v <- rlogis(d*(d-1)/2)
  L <- cvine_cholesky(v, d, eta = 1)
  R <- L %*% t(L)
  expect_equal(diag(R), rep(1, d))
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})

test_that("rlkj_cvine returns correlation matrices", {
  R <- rlkj_cvine(1, 3, eta = 2)
  expect_equal(dim(R), c(3,3))
  expect_equal(diag(R), rep(1,3))
  expect_true(all(eigen(R, only.values = TRUE)$values > 0))
})
