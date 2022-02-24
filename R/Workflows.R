#' Workflow for Random Forests
#'
#' @param form A formula object.
#' @param train A train dataset. Data.frame or tibble object.
#' @param test A test dataset. Data.frame or tibble object.
#' @param ... Additional paramters which can be passed into internal model.
#'
#' @return A list containing true values and predictions.
#' @export
#'
#' @examples
wf.ranger <- function(form, train, test, ...){

  m <- ranger(formula = form, data =  train, ...)
  preds <- predict(m, test)$predictions

  r <- list(trues = responseValues(form, test),
            preds = preds)
  return(r)
}

#' Workflow for xgboost
#'
#' @param form A formula object.
#' @param train A train dataset. Data.frame or tibble object.
#' @param test A test dataset. Data.frame or tibble object.
#' @param ... Additional paramters which can be passed into internal model.
#'
#' @return A list containing true values and predictions.
#' @export
#'
#' @examples

wf.xgboost <- function(formula, train, test,...){

  t <- formula[[2]]

  xgb_train_m <- data.matrix(dplyr::select(train, -t))
  xgb_test_m <- data.matrix(dplyr::select(test, -t))

  label_train <- pull(train, all_of(t))
  label_test <- pull(test, all_of(t))

  xgb_train <- xgb.DMatrix(data = xgb_train_m, label = label_train)
  xgb_test <- xgb.DMatrix(data = xgb_test_m, label = label_test)



  m <- xgboost(data = xgb_train,
               verbose = 0,
               booster = "dart",
               ...
  )

  preds <- predict(m, xgb_test)

  r <- list(trues = as.vector(responseValues(formula, test)),
            preds = preds
  )
  return(r)

}

#' Workflow for xgboost optimised with SERA
#'
#' @param form A formula object.
#' @param train A train dataset. Data.frame or tibble object.
#' @param test A test dataset. Data.frame or tibble object.
#' @param ... Additional parameters which can be passed into internal model.
#'
#' @return A list containing true values and predictions.
#' @export
#'
#' @examples

wf.xSERAgboost <- function(formula, train, test,...){

  t <- formula[[2]]

  xgb_train_m <- data.matrix(dplyr::select(train, -all_of(t)))
  xgb_test_m <- data.matrix(dplyr::select(test, -all_of(t)))

  label_train <- pull(train, all_of(t))
  label_test <- pull(test, all_of(t))

  xgb_train <- xgboost::xgb.DMatrix(data = xgb_train_m, label = label_train)
  xgb_test <- xgboost::xgb.DMatrix(data = xgb_test_m, label = label_test)


  m <- xgboost::xgboost(data = xgb_train,
                        verbose = 0,
                        booster = "dart",
                        objective = xgboostsera,
                        #eval_metric = seraerr,
                        ...
  )

  preds <- predict(m, xgb_test)

  r <- list(trues = as.vector(responseValues(formula, test)),
            preds = preds
  )
  return(r)

}

#' Workflow for Gradient Boost optimised with SERA.
#'
#' @param form A formula object.
#' @param train A train dataset. Data.frame or tibble object.
#' @param test A test dataset. Data.frame or tibble object.
#' @param ... Additional paramters which can be passed into internal model.
#'
#' @return A list containing true values and predictions.
#' @export
#'
#' @examples

wf.SERAGradientBoost <- function(form, train, test,...){

  preds <- SERAGradientBoost(formula = form, train, test, ...)
  r <- list(trues = responseValues(form, test),
            preds = preds)
  return(r)
}


#' Workflow for Gradient Tree Boost optimised with SERA.
#'
#' @param form A formula object.
#' @param train A train dataset. Data.frame or tibble object.
#' @param test A test dataset. Data.frame or tibble object.
#' @param ... Additional paramters which can be passed into internal model.
#'
#' @return A list containing true values and predictions.
#' @export
#'
#' @examples

wf.SERAGradientTreeBoost <- function(form, train, test,...){

  preds <- SERAGradientTreeBoost(formula = form, train, test, ...)

  r <- list(trues = responseValues(form, test),
            preds = preds)
  return(r)
}

#' Workflow for Gradient Boost
#'
#' @param form A formula object.
#' @param train A train dataset. Data.frame or tibble object.
#' @param test A test dataset. Data.frame or tibble object.
#' @param ... Additional paramters which can be passed into internal model.
#'
#' @return A list containing true values and predictions.
#' @export
#'
#' @examples

wf.GradientBoost <- function(form, train, test,...){

  preds <- GradientBoost(formula = form, train, test, ...)

  r <- list(trues = responseValues(form, test),
            preds = preds)
  return(r)
}

#' Workflow for Gradient Tree Boost
#'
#' @param form A formula object.
#' @param train A train dataset. Data.frame or tibble object.
#' @param test A test dataset. Data.frame or tibble object.
#' @param ... Additional paramters which can be passed into internal model.
#'
#' @return A list containing true values and predictions.
#' @export
#'
#' @examples


wf.GradientTreeBoost <- function(form, train, test,...){

  preds <- GradientTreeBoost(formula = form, train, test, ...)

  r <- list(trues = responseValues(form, test),
            preds = preds)
  return(r)
}

#' Workflow for Light Gradient Boost.
#'
#' @param form A formula object.
#' @param train A train dataset. Data.frame or tibble object.
#' @param test A test dataset. Data.frame or tibble object.
#' @param ... Additional paramters which can be passed into internal model.
#'
#' @return A list containing true values and predictions.
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

  label_train <- pull(train, all_of(t))

  lgbm_train <- lgb.Dataset(data = lgbm_train_m, label = label_train)

  pars <- list(...)

  m <- lightgbm::lightgbm(data = lgbm_train,
                          label = label_train,
                          obj = "regression",
                          verbose = -1,
                          boosting = "dart",
                          params = pars,
                          min_data = min_data)



  preds <- predict(m, lgbm_test_m)

  r <- list(trues = as.vector(responseValues(formula, test)),
            preds = preds)
  return(r)

}

#' Workflow for Light Gradient Boost optimised with SERA.
#'
#' @param form A formula object.
#' @param train A train dataset. Data.frame or tibble object.
#' @param test A test dataset. Data.frame or tibble object.
#' @param ... Additional paramters which can be passed into internal model.
#'
#' @return A list containing true values and predictions.
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

  label_train <- pull(train, all_of(t))

  lgbm_train <- lgb.Dataset(data = lgbm_train_m, label = label_train)

  pars <- list(...)

  m <- lightgbm::lightgbm(data = lgbm_train,
                          label = label_train,
                          obj = lgbmsera,
                          verbose = -1,
                          boosting = "dart",
                          params = pars,
                          min_data = min_data)



  preds <- predict(m, lgbm_test_m)

  r <- list(trues = as.vector(responseValues(formula, test)),
            preds = preds)
  return(r)

}
