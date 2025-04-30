#' lkj
#' Return the Cholesky factor of the correlation matrix from random parameter
#'
#' @param k integer. The dimension of the matrix
#' @param y The number of parameters to map.
#' The number of parameters should be is k*(k-1)/2
#'
#' @returns A Square correlation matrix
#' @export
#'
#' @examples
#' R_matrix(2, 0.1)
L_matrix <- function(k, y){
  t(.lkj_cmit(k = k, y))
  # t(w) %*% w
}
