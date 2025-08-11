#' lkj
#' Return the parameters  of a Cholesky factor of the correlation matrix
#'
#' @param L Matrix. The Cholesky factor matrix
#' The number of parameters should be is k*(k-1)/2
#'
#' @returns A vector of parameters of k*(k-1)/2 lenght
#' @export
#'
#' @examples
#' L <- L_matrix(2, 0.1)
#' L_par(L)
L_par <- function(L){
  .lkj_par( t(L) )
}
