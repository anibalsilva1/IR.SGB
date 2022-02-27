#' Returns the workflow for each model which has the higher posterior probability of being
#' practically better than another workflow.
#'
#' @param bayesresults A data.frame returned from \code{\link{get_all_bayes_models_results.R}}
#'
#' @return A character vector with the best workflow for all the considered models.
#' @export
#'
#' @examples

get_best_oracles_for_each_model <- function(bayesresults){

  models <- names(bayesresults)
  best <- sapply(models,
                 FUN = function(model) bayesresults[[model]] %>% group_by(oracle) %>%
                   summarise(prob = mean(oracleProb), .groups = "drop") %>%
                   slice_max(prob, n = 1, with_ties = FALSE) %>%
                   pull(oracle)
  )
  return(best)
}
