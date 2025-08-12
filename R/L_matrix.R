#' lkj
#' Return the Cholesky factor of the correlation matrix from random parameter
#'
#' @param y The number of parameters to map.
#' The number of parameters should be is k*(k-1)/2
#'
#' @returns A Square correlation matrix
#' @export
#'
#' @examples
#' L_matrix(0.1)
L_matrix <- function(y){
  #inverse of k = n*(n-1)/2
  f <- function(n) 0.5 * (1 + sqrt(8*n + 1))
  k <- f(length(y))
  t(.lkj_cmit(k = k, y))
  # t(w) %*% w
}
