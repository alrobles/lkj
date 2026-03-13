#' Deterministic mapping from unconstrained reals to a Cholesky factor or
#' correlation matrix (LKJ distribution)
#'
#' Maps a vector of unconstrained real numbers to the Cholesky factor \eqn{L}
#' (or the correlation matrix \eqn{R = L L^\top}) of a correlation matrix whose
#' density under the LKJ distribution is proportional to
#' \eqn{\det(R)^{\eta - 1}}.
#'
#' Two construction methods are supported:
#' \describe{
#'   \item{\code{"cvine"}}{C-vine method (Lewandowski et al. 2009, Sec. 2.4).
#'     Requires \eqn{d(d-1)/2} parameters.}
#'   \item{\code{"onion"}}{Extended onion method (Ghosh & Henderson 2003;
#'     Lewandowski et al. 2009). Requires \eqn{d(d+1)/2 - 2} parameters.}
#' }
#'
#' @param v Numeric vector of unconstrained real parameters. The required length
#'   depends on \code{d} and \code{method} (see Details).
#' @param d Positive integer. Dimension of the target correlation matrix.
#' @param eta Positive numeric. LKJ shape parameter. \code{eta = 1} (default)
#'   gives the uniform distribution over correlation matrices; larger values
#'   concentrate mass near the identity.
#' @param method Character string: \code{"cvine"} (default) or \code{"onion"}.
#' @param output Character string: \code{"L"} to return the lower-triangular
#'   Cholesky factor, or \code{"R"} (default) to return the full correlation
#'   matrix \eqn{R = L L^\top}.
#'
#' @return A \code{d x d} matrix: the Cholesky factor \eqn{L} if
#'   \code{output = "L"}, or the correlation matrix \eqn{R} if
#'   \code{output = "R"}.
#'
#' @details
#' Required parameter lengths:
#' \itemize{
#'   \item \code{"cvine"}: \eqn{d(d-1)/2}
#'   \item \code{"onion"}: \eqn{d(d+1)/2 - 2} for \eqn{d \geq 2}; 0 for
#'     \eqn{d = 1}.
#' }
#'
#' @references
#' Lewandowski, D., Kurowicka, D., & Joe, H. (2009). Generating random
#' correlation matrices based on vines and extended onion method.
#' \emph{Journal of Multivariate Analysis}, 100(9), 1989–2001.
#'
#' @examples
#' # C-vine: d=4 needs 4*(4-1)/2 = 6 parameters
#' v <- c(-0.5, 0.3, 1.2, -0.8, 0.1, 0.6)
#' R <- lkj(v, d = 4, eta = 1, method = "cvine")
#' print(round(R, 4))
#'
#' # Onion: d=4 needs 4*(4+1)/2 - 2 = 8 parameters
#' v2 <- rep(0, 8)
#' R2 <- lkj(v2, d = 4, eta = 1, method = "onion")
#' print(round(R2, 4))
#'
#' # Return Cholesky factor instead of full matrix
#' L <- lkj(v, d = 4, eta = 1, method = "cvine", output = "L")
#'
#' @export
lkj <- function(v, d, eta = 1,
                method = c("cvine", "onion"),
                output = c("R", "L")) {
  method <- match.arg(method)
  output <- match.arg(output)

  L <- if (method == "cvine") {
    cvine_cholesky(v, d, eta)
  } else {
    onion_cholesky(v, d, eta)
  }

  if (output == "L") L else L %*% t(L)
}
