#' Percentage difference
#'
#' @description Evaluates the percentage difference between models that were
#' optimised with standard metrics vs models that were optimised with SERA.
#'
#' @param preds A \code{list} with predictions for each dataset. Each element of
#' the \code{list} should be a \code{data.frame} where each column is a vector with
#' the predictions of a given model.
#' @param models A character \code{vector} with the name of the models.
#'
#' @return A \code{data.frame} with the percentual differences for each model.
#' @export
#'
#' @examples
percentual_diff <- function(preds, models){

  n <- length(preds)
  step <-  0.001
  s <- seq(0, 1, step)
  m <- length(s)

  res <- data.frame(phi = s)

  for(model in models){

    print(paste("Iteration over model ", model))

    errs_std <- c()
    errs_nonstd <- c()
    errs_ds <- matrix(nrow = m, ncol = n)
    errs <- c()

    for(i in 1:n){

      ds <- preds[[n]]$preds
      ph.ctrl <- preds[[n]]$p.ctrl

      tr <- ds$trues

      m_preds <- ds %>% dplyr::select(model) %>% pull()
      sm_preds <- ds %>% dplyr::select(matches(str_c(model, "_SERA"))) %>% pull()

      errs_std <- sera(tr, m_preds, ph = ph.ctrl, return.err = T)$errors
      errs_nonstd <- sera(tr, sm_preds, ph = ph.ctrl, return.err = T)$errors
      errs_ds[, i] <- 100 * (errs_nonstd - errs_std) / errs_std
    }

    errs <- rowMeans(errs_ds)

    res[model] <- errs

  }

  return(res)
}
