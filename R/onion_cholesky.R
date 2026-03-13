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

  onion_cholesky_cpp(v, d, eta)
}
