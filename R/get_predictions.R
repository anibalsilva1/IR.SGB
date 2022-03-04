#' Performs predictions for the best models given a list of datasets.
#'
#' @param datasets A list containing the datasets with their information.
#' @param bestmodels A dataframe with the best models.
#'
#' @return A list with the predictions for each dataset. This predictions
#' will be inside a dataframe, where each column will be the predictions of each model.
#' @export
#' @import doParallel
#' @importFrom foreach foreach %dopar%
#' @examples
get_predictions <- function(datasets, bestmodels){

  nD <- names(datasets)
  nWF <- pull(distinct(bestmodels, workflow))

  preds_data <- list()

  for(i in 1:length(nD)){

    form <- datasets[[i]]$formula
    train <- datasets[[i]]$train
    test <- datasets[[i]]$test
    target <- form[[2]]


    p.ctrl = phi.control(pull(train,target))

    df <- tibble(trues = pull(test,target))

    print(paste0("Gettings results for ", nD[i], " data"))

    preds_wf <- foreach(wf = 1:length(nWF),
                        .final = function(wf) setNames(wf, nWF)
    ) %dopar% {

      if(nWF[wf] != "wf.ranger")
        params <- bestmodels[bestmodels$dataset == nD[i] &
                               bestmodels$workflow == nWF[wf], ]$params[[1]]
      else
        params <- bestmodels[bestmodels$dataset == nD[i] &
                               bestmodels$workflow == nWF[wf], ]$params[[1]]$learner.pars

      preds <- do.call(what = nWF[wf],
                       args = c(list(form, train, test), params))$preds
    }


    preds_df <- as_tibble(do.call("cbind", preds_wf))
    df <- bind_cols(df, preds_df)

    preds_data[[nD[i]]] <- list(preds = df,
                                p.ctrl = p.ctrl)

  }
  return(preds_data)
}
