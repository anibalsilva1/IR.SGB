#' Get times
#'
#' @description Evaluates execution time for a set of models in a given set of datasets.
#'
#' @param datasets A \code{list} of datasets containing a \code{formula} object,
#' a \code{data.frame} object with a train set and a \code{data.frame} object with
#' a test set.
#' @param parameters A \code{data.frame} object with a set of models and parameters.
#'
#' @return A \code{list} where each element is a \code{data.frame} with the respective
#' execution times.
#' @export
#'
#' @examples
get_times <- function(datasets, parameters){

  nD <- names(datasets)
  nWF <- pull(distinct(parameters, workflow))
  preds_data <- list()

  for(i in 1:length(nD)){

    form <- datasets[[i]]$formula
    train <- datasets[[i]]$train
    test <- datasets[[i]]$test
    target <- form[[2]]


    p.ctrl = phi.control(pull(train,target))

    df <- tibble()

    print(paste0("Gettings results for ", nD[i], " data"))

    preds_wf <- foreach(wf = 1:length(nWF),
                        .final = function(wf) setNames(wf, nWF)
    ) %dopar% {

      if(nWF[wf] != "wf.ranger")
        params <- parameters[parameters$dataset == nD[i] &
                               parameters$workflow == nWF[wf], ]$params[[1]]
      else
        params <- parameters[parameters$dataset == nD[i] &
                               parameters$workflow == nWF[wf], ]$params[[1]]$learner.pars

      preds <- do.call(what = nWF[wf],
                       args = c(list(form, train, test, time = T), params))$time
    }

    mg <- do.call("rbind", preds_wf)
    df <- as_tibble(mg)
    print(df)

    df <- df %>% mutate(model = nWF) %>% relocate(model)
    colnames(df) <- c("model","train_time", "test_time")

    preds_data[[nD[i]]] <- df

  }
  return(preds_data)
}
