#' Determines predictions of a given dataset using the Gradient Boost for Regression using
#' SERA as an optimisation loss function.
#'
#'
#' @param formula A formula object.
#' @param train The training dataset. An data.frame or tibble object.
#' @param test The test dataset. An data.frame or tibble object.
#' @param maxIter The maximum number of iterations.
#' @param eta Learning rate.
#' @param lambda Ridge regression parameter.
#' @param maxdepth Max depth of a tree.
#' @param verbose Prints out the error across iterations (if 1).
#'
#' @return A numeric vector with predictions.
#' @export
#'
#' @examples
#'


library(IRon)
library(rpart)
library(treeClust)

SERAGradBoost <- function(formula,
                          train,
                          test,
                          maxIter = 100,
                          eta = 0.01,
                          lambda = 10,
                          maxdepth = 5,
                          verbose = 0){

  ### FIX ROWNAMES BECAUSE OF UNORDERED LEAF NAMES
  rownames(train) <- 1:nrow(train)
  rownames(test) <- 1:nrow(test)

  target <- formula[[2]]

  y <- dplyr::pull(train, target)
  X <- dplyr::select(train, -target)

  stumps <- list()
  gammas <- c()


  phi.trues <- phi(y)
  ser_m <- sera_min(y, phi.trues)

  s <- 0.001
  step <- seq(0,1,s)
  N <- length(step)

  sigmas <- get_sigma(phi.trues, step)
  sigmas <- sigmas/N


  rows = nrow(train)

  weakpreds   <- c()
  strongpreds <- c()
  pseudo_res  <- c()

  erra <- c()

  t = 1
  F_0 <- ser_m
  strongpreds <- rep(F_0, rows)

  while (t <= maxIter) {

    X$pseudo_res <- 2*sigmas*(y - strongpreds)

    resWeak <- rpart(formula = pseudo_res ~ .,
                     data = X,
                     maxdepth = maxdepth)

    weakpreds <- predict(resWeak, X)
    stumps[[t]] <- resWeak

    errs_num <- sapply(step, FUN = function(i) sum(weakpreds[phi.trues >= i]*(y[phi.trues >= i] - strongpreds[phi.trues >= i])))
    errs_den <- sapply(step, FUN = function(i) sum(weakpreds[phi.trues >= i]^2) + lambda)

    areas_num <- sum(sapply(2:length(step), FUN = function(x) s * (errs_num[x-1] + errs_num[x])/ 2 ))
    areas_den <- sum(sapply(2:length(step), FUN = function(x) s * (errs_den[x-1] + errs_den[x])/ 2 ))

    gammas[t] <- areas_num/areas_den

    strongpreds <- strongpreds + eta * gammas[t] * weakpreds
    X <- dplyr::select(X, -pseudo_res)



    if(verbose == 1){
      erra[t] <- sera(trues = y, preds = strongpreds, phi.trues = phi.trues)
      print(paste0("Iteration: ", t, " SERA: ", erra[t]))
    }

    t = t+1

  }

  m = nrow(test)
  n = length(stumps)
  preds <- matrix(0, nrow = m, ncol = n)

  for (i in 1:n) {

    preds[, i] <- predict(stumps[[i]], test)

  }
  preds <- F_0 + sapply(1:m, FUN = function(i) eta * sum(gammas * preds[i ,]))

  return(preds)
}
