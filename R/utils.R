#' Custom adaptation of SERA to xgboost.
#'
#' @param preds A numeric vector of predictions.
#' @param dtrain A xgb.DMatrix object.
#'
#' @return Returns a list containing the gradient and the hessian of SERA.
#' @export
#'
#' @examples
xgboostsera <- function(preds, dtrain){


  s <- 0.001
  labels <- xgboost::getinfo(dtrain, "label")

  step <- seq(0,1,s)
  N <- length(step)

  phi.ctrl <- phi.control(labels)
  phi.trues <- phi(labels, phi.ctrl)
  sigmas <- get_sigma(phi.trues, step)
  sigmas <- sigmas/N


  grad <- 2*sigmas*(preds - labels)
  hess <- 2*sigmas

  return(list(grad = grad, hess = hess))
}

#' Custom adaptation of SERA to Light GBM.
#'
#' @param preds A numeric vector of predictions.
#' @param dtrain A lgb.Dataset object.
#'
#' @return  Returns a list containing the gradient and the hessian of SERA.
#' @export
#'
#' @examples
lgbmsera <- function(preds, dtrain){

  labels <- lightgbm::get_field(dtrain, "label")


  s <- 0.001
  step <- seq(0,1,s)
  N <- length(step)

  #
  phi.ctrl <- phi.control(labels)
  phi.trues <- phi(labels, phi.ctrl)
  sigmas <- get_sigma(phi.trues, step)
  sigmas <- sigmas/N

  #grad <- sera_deriv(labels, preds, phi = phi.trues)
  grad <- 2*sigmas*(preds - labels)
  hess <- 2*sigmas

  return(list(grad = grad, hess = hess))
}
