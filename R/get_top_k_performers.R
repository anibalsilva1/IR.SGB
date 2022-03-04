#' Given a data.frame with all the workflows for each data.frame, returns a character
#' vector with \code{k} workflows that had the lowest error across all datasets.
#'
#' @param res A data.frame returned from \code{\link{get_avg_scores_by_metric}}.
#' @param k A numeric with the number of workflows to return.
#' @param m.name A character with the name of the model.
#' @param metric A character with the metric.
#'
#' @return A character vector with the top-k workflows.
#' @export
#'
#' @examples
get_top_k_performers <- function(res, k = 3, m.name, metric = "mse"){

  res <- res %>% tidyr::separate(workflow,
                          into = c("model", "version"),
                          sep = ".v") %>%
    dplyr::filter(model == m.name)

    topMSE <- res %>% group_by(dataset) %>%
      slice_min(avg_mse, n = k) %>%
      ungroup() %>%
      count(version) %>%
      slice_max(version, n = k) %>%
      select(-n)

    topSERA <- res %>% group_by(dataset) %>%
      slice_min(avg_sera, n = k) %>%
      ungroup() %>%
      count(version) %>%
      slice_max(version, n = k) %>%
      select(-n)

    top <- bind_rows(topMSE, topSERA) %>% distinct(version) %>% pull()

    res <- sapply(top, FUN = function(i) paste0(m.name, ".v", i))

  return(res)
}
