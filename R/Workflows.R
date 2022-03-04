#' XGBoost workflow
#'
#' @description Workflow for XGBoost.
#'
#' @param formula A \code{formula} object.
#' @param train \code{data.frame} or \code{tibble} object with the training set.
#' @param test \code{data.frame} or \code{tibble} object with the test set.
#' @param ... Additional parameters which can be passed into internal model.
#'
#' @return A \code{list} containing true values and predictions.
#' @export
#'
#' @examples

wf.xgboost <- function(formula, train, test,...){

  t <- formula[[2]]

  xgb_train_m <- data.matrix(dplyr::select(train, -t))
  xgb_test_m <- data.matrix(dplyr::select(test, -t))

  label_train <- dplyr::pull(train, all_of(t))
  label_test <- dplyr::pull(test, all_of(t))

  xgb_train <- xgboost::xgb.DMatrix(data = xgb_train_m, label = label_train)
  xgb_test <-  xgboost::xgb.DMatrix(data = xgb_test_m, label = label_test)



  m <- xgboost::xgboost(data = xgb_train,
                        verbose = 0,
                        booster = "dart",
                        ...
  )

  preds <- predict(m, xgb_test)

  r <- list(trues = as.vector(performanceEstimation::responseValues(formula, test)),
            preds = preds
  )
  return(r)

}

#' XGBoost(SERA) workflow
#'
#' @description Workflow for XGBoost optimised with SERA.
#'
#' @param formula A \code{formula} object.
#' @param train \code{data.frame} or \code{tibble} object with the training set.
#' @param test \code{data.frame} or \code{tibble} object with the test set.
#' @param ... Additional parameters which can be passed into internal model.
#'
#' @return A \code{list} containing true values and predictions.
#' @export
#'
#' @examples

wf.xSERAgboost <- function(formula, train, test,...){

  t <- formula[[2]]

  xgb_train_m <- data.matrix(dplyr::select(train, -all_of(t)))
  xgb_test_m <- data.matrix(dplyr::select(test, -all_of(t)))

  label_train <- dplyr::pull(train, all_of(t))
  label_test <-  dplyr::pull(test, all_of(t))

  xgb_train <- xgboost::xgb.DMatrix(data = xgb_train_m, label = label_train)
  xgb_test <- xgboost::xgb.DMatrix(data = xgb_test_m, label = label_test)


  m <- xgboost::xgboost(data = xgb_train,
                        verbose = 0,
                        booster = "dart",
                        objective = xgboostsera,
                        ...
  )

  preds <- predict(m, xgb_test)

  r <- list(trues = as.vector(performanceEstimation::responseValues(formula, test)),
            preds = preds
  )
  return(r)

}

#' GB(SERA) workflow
#'
#' @description Workflow for Gradient Boosting optimised with SERA.
#'
#' @param formula A \code{formula} object.
#' @param train \code{data.frame} or \code{tibble} object with the training set.
#' @param test \code{data.frame} or \code{tibble} object with the test set.
#' @param ... Additional parameters which can be passed into internal model.
#'
#' @return A \code{list} containing true values and predictions.
#' @export
#'
#' @examples

wf.SERAGradientBoost <- function(formula, train, test,...){

  preds <- SERAGradientBoost(formula = formula, train, test, ...)
  r <- list(trues = performanceEstimation::responseValues(formula, test),
            preds = preds)
  return(r)
}


#' GBRT (SERA)
#'
#' @description Workflow for Gradient Boosting Regression Trees optimised with SERA.
#'
#' @param formula A \code{formula} object.
#' @param train \code{data.frame} or \code{tibble} object with the training set.
#' @param test \code{data.frame} or \code{tibble} object with the test set.
#' @param ... Additional parameters which can be passed into internal model.
#'
#' @return A \code{list} containing true values and predictions.
#' @export
#'
#' @examples

wf.SERAGradientTreeBoost <- function(formula, train, test,...){

  preds <- SERAGradientTreeBoost(formula = formula, train, test, ...)

  r <- list(trues = performanceEstimation::responseValues(formula, test),
            preds = preds)
  return(r)
}

#' GB
#'
#' @description Workflow for Gradient Boosting.
#'
#' @param formula A \code{formula} object.
#' @param train \code{data.frame} or \code{tibble} object with the training set.
#' @param test \code{data.frame} or \code{tibble} object with the test set.
#' @param ... Additional parameters which can be passed into internal model.
#'
#' @return A \code{list} containing true values and predictions.
#' @export
#'
#' @examples

wf.GradientBoost <- function(formula, train, test,...){

  preds <- GradientBoost(formula = formula, train, test, ...)

  r <- list(trues = performanceEstimation::responseValues(formula, test),
            preds = preds)
  return(r)
}

#' GBRT
#'
#' @description Workflow for Gradient Boosting Regression Trees.
#'
#' @param formula A \code{formula} object.
#' @param train \code{data.frame} or \code{tibble} object with the training set.
#' @param test \code{data.frame} or \code{tibble} object with the test set.
#' @param ... Additional parameters which can be passed into internal model.
#'
#' @return A \code{list} containing true values and predictions.
#' @export
#'
#' @examples


wf.GradientTreeBoost <- function(formula, train, test,...){

  preds <- GradientTreeBoost(formula = formula, train, test, ...)

  r <- list(trues = performanceEstimation::responseValues(formula, test),
            preds = preds)
  return(r)
}

#' LGB
#'
#' @description Workflow for Light Gradient Boosting Machines.
#'
#' @param formula A \code{formula} object.
#' @param train \code{data.frame} or \code{tibble} object with the training set.
#' @param test \code{data.frame} or \code{tibble} object with the test set.
#' @param ... Additional parameters which can be passed into internal model.
#'
#' @return A \code{list} containing true values and predictions.
#' @export
#'
#' @examples

wf.LGBM <- function(formula, train, test, ...){

  t <- formula[[2]]

  if(t == "C_peptide")
    min_data = 15
  else
    min_data = 20

  lgbm_train_m <- data.matrix(dplyr::select(train, -all_of(t)))
  lgbm_test_m <- data.matrix(dplyr::select(test, -all_of(t)))

  label_train <- dplyr::pull(train, all_of(t))

  lgbm_train <- lightgbm::lgb.Dataset(data = lgbm_train_m, label = label_train)

  pars <- list(...)

  m <- lightgbm::lightgbm(data = lgbm_train,
                          label = label_train,
                          obj = "regression",
                          verbose = -1,
                          boosting = "dart",
                          params = pars,
                          min_data = min_data)



  preds <- predict(m, lgbm_test_m)

  r <- list(trues = as.vector(performanceEstimation::responseValues(formula, test)),
            preds = preds)
  return(r)

}

#' LGB(SERA)
#'
#' @description Workflow for Light Gradient Boosting optimised with SERA.
#'
#' @param formula A \code{formula} object.
#' @param train \code{data.frame} or \code{tibble} object with the training set.
#' @param test \code{data.frame} or \code{tibble} object with the test set.
#' @param ... Additional parameters which can be passed into internal model.
#'
#' @return A \code{list} containing true values and predictions.
#' @export
#'
#' @examples

wf.LGBMSERA <- function(formula, train, test, ...){

  t <- formula[[2]]

  if(t == "C_peptide")
    min_data = 15
  else
    min_data = 20

  lgbm_train_m <- data.matrix(dplyr::select(train, -all_of(t)))
  lgbm_test_m <- data.matrix(dplyr::select(test, -all_of(t)))

  label_train <- dplyr::pull(train, all_of(t))

  lgbm_train <- lightgbm::lgb.Dataset(data = lgbm_train_m, label = label_train)

  pars <- list(...)

  m <- lightgbm::lightgbm(data = lgbm_train,
                          label = label_train,
                          obj = lgbmsera,
                          verbose = -1,
                          boosting = "dart",
                          params = pars,
                          min_data = min_data)



  preds <- predict(m, lgbm_test_m)

  r <- list(trues = as.vector(performanceEstimation::responseValues(formula, test)),
            preds = preds)
  return(r)

}
