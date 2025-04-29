#' lkj
#' Return correlation matrix from random parameter
#'
#' @param k integer. The dimension of the matrix
#' @param y The number of parameters to map.
#' The number of parameters should be is k*(k-1)/2
#'
#' @returns A Square correlation matrix
#' @export
#'
#' @examples
#' L(2, 1)
L <- function(k, y){
  w <- .lkj_cmit(k = k, y)
  t(w) %*% w
}
