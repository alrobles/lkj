#' Generate random correlation matrices from the LKJ distribution
#'
#' Draws \code{n} random correlation matrices of dimension \code{d} from the
#' LKJ distribution, whose density is proportional to \eqn{\det(R)^{\eta - 1}}.
#' Dispatches to either the C-vine or extended onion construction method.
#'
#' @param n Positive integer. Number of matrices to generate.
#' @param d Positive integer. Dimension of the matrices.
#' @param eta Positive numeric. LKJ shape parameter. \code{eta = 1} (default)
#'   gives the uniform distribution over correlation matrices; larger values
#'   concentrate mass near the identity.
#' @param method Character string: \code{"cvine"} (default) or \code{"onion"}.
#' @param output Character string: \code{"R"} (default) to return correlation
#'   matrices, or \code{"L"} to return the lower-triangular Cholesky factors.
#'
#' @return If \code{n = 1}, a \code{d x d} matrix.
#'   If \code{n > 1}, an array of dimension \code{c(d, d, n)}.
#'
#' @references
#' Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random
#' correlation matrices based on vines and extended onion method.
#' \emph{Journal of Multivariate Analysis}, 100(9), 1989–2001.
#'
#' @examples
#' # Single 3x3 correlation matrix via C-vine
#' R <- rlkj(1, d = 3, eta = 2)
#'
#' # Ten 4x4 matrices via onion method
#' arr <- rlkj(10, d = 4, eta = 1, method = "onion")
#'
#' # Return Cholesky factors instead of correlation matrices
#' L_arr <- rlkj(5, d = 3, eta = 1, method = "cvine", output = "L")
#'
#' @export
rlkj <- function(n, d, eta = 1,
                 method = c("cvine", "onion"),
                 output = c("R", "L")) {
  method <- match.arg(method)
  output <- match.arg(output)

  stopifnot(
    "n must be a positive integer" = n == as.integer(n) && n >= 1L,
    "d must be a positive integer" = d == as.integer(d) && d >= 1L,
    "eta must be a positive number" = is.numeric(eta) && length(eta) == 1L && eta > 0
  )

  draw_one <- function() {
    if (method == "cvine") {
      m <- d * (d - 1L) / 2L
      v <- stats::rlogis(m)
      cvine_cholesky(v, d, eta)
    } else {
      m <- if (d <= 1L) 0L else if (d == 2L) 1L else d * (d + 1L) / 2L - 2L
      v <- stats::rlogis(m)
      onion_cholesky(v, d, eta)
    }
  }

  if (n == 1L) {
    L <- draw_one()
    return(if (output == "L") L else L %*% t(L))
  }

  out <- array(0, dim = c(d, d, n))
  for (i in seq_len(n)) {
    L <- draw_one()
    out[, , i] <- if (output == "L") L else L %*% t(L)
  }
  out
}
