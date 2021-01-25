#' Wilcoxon Sign Rank Test Statistic Exact Distribution
#'
#' @name W_stat
#' @description This function allows the user to find the probability
#' values from the exact distribution of W, Bickel and Doksum(1973).
#' The exact P(W=x), P(W<=x), P(W>=x) values is found via an exhaustive enumeration
#' of the possible permutations of data with size n.
#'
#' @usage W_stat(n , test_stat, side = c('geq','leq','eq'))
#' @param n Size of data or Number of observations
#' @param test_stat The x value specified in  P(W=x), P(W<=x), P(W>=x)
#' @param side The tails of exact probability the user wants to compute  e.g.
#' 'eq' = P(W=x), 'leq' = P(W<=x), 'geq' = 'P(W>=x)
#' @return The exact probability values as specified.
#' @examples
#' W_stat(n=5, test_stat = 3, side = 'leq')
#' @export
#'
#'

Sys.setenv('_R_CHECK_SYSTEM_CLOCK_' = 0)

W_stat = function(n, test_stat , side = c('geq','leq','eq')){
  mat = expand.grid(rep(list(c(-1,1)),(n)))
  names(mat) = 1:n ; vec = 1:n
  mat = sweep(mat, MARGIN=2, vec, `*`)
  positive = 1*(mat>0)
  mat= cbind(mat, "positive sum" =apply(mat*positive, 1, sum))
  if(side == 'geq'){
    num_combi = sum(mat$`positive sum` >= test_stat)
  }
  else if(side == 'leq'){
    num_combi = sum(mat$`positive sum` <= test_stat)
  }
  else if(side == 'eq'){num_combi = sum(mat$`positive sum` == test_stat)}
  total = nrow(mat)
  prob = num_combi/total
  return(prob)
}

