#' Calculates the posterior probabilities of a model A (left) being practically better then a model B (right).
#'
#' @param diffVector A numeric vector with error differences.
#' @param rope_min A numeric value with the minimum of ROPE.
#' @param rope_max A numeric value with the maximum of ROPE.
#'
#' @return A list with the posterior probabilities.
#' @export
#'
#' @examples

BayesianSignTest <- function(diffVector,rope_min,rope_max) {

  samples <- 3000

  #build the vector 0.5 1 1 ....... 1
  weights <- c(0.5,rep(1,length(diffVector)))

  #add the fake first observation in 0
  diffVector <- c (0, diffVector)


  #for the moment we implement the sign test. Signedrank will follows
  probLeft <- mean (diffVector < rope_min)
  probRope <- mean (diffVector > rope_min & diffVector < rope_max)
  probRight <- mean (diffVector > rope_max)

  alpha <- c(probLeft,probRope,probRight)

  alpha <- alpha+0.0001

  alpha.res <- colMeans(MCMCpack::rdirichlet(30000, alpha))

  results = list("probLeft"=alpha.res[1], "probRope"=alpha.res[2],
                 "probRight"=alpha.res[3])

  return (results)

}
