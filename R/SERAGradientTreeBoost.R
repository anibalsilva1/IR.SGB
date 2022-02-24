#' Determines predictions of a given dataset using the Gradient Tree Boost for Regression using
#' SERA as an optimisation loss function.
#'
#' @param formula A formula object.
#' @param train The training dataset. An data.frame or tibble object.
#' @param test The test dataset. An data.frame or tibble object.
#' @param maxIter The maximum number of iterations.
#' @param eta Learning rate.
#' @param verbose Prints out the error across iterations (if 1)
#'
#' @return A numeric vector with predictions.
#' @export
#'
#' @examples
#'

library(IRon)
library(rpart)
library(treeClust)

SERAGradientTreeBoost <- function(formula,
                                  train,
                                  test,
                                  maxIter = 100,
                                  eta = 0.01,
                                  verbose = 0){

  rownames(train) <- 1:nrow(train)
  rownames(test) <- 1:nrow(test)

  target <- formula[[2]]

  y <- dplyr::pull(train, target)
  X <- dplyr::select(train, -target)
  df <- X

  phi.ctrl <- phi.control(y)

  phi.trues <- phi(y, phi.ctrl)

  s <- 0.001
  step <- seq(0,1,s)
  N <- length(step)

  sigmas <- get_sigma(phi.trues, step)/N

  rows = nrow(train)

  stumps <- list()
  gammas <- list()
  leaves <- list()

  weakpreds   <- c()
  strongpreds <- c()
  pseudo_res  <- c()
  erra <- c()

  t = 1
  F_0 = sera_min(y, phi.trues)
  strongpreds <- rep(F_0, rows)


  while (t <= maxIter) {

    df$pseudo_res <- 2*sigmas*(y - strongpreds)

    resWeak <- rpart(formula = pseudo_res ~ ., data = df)
    stumps[[t]] <- resWeak
    leaves[[t]] <- resWeak$where

    ers_num <- sapply(unique(leaves[[t]]), FUN = function(leaf)
      sapply(step, FUN = function(i) sum(y[phi.trues >= i & leaves[[t]] == leaf] - strongpreds[phi.trues >= i & leaves[[t]] == leaf])))

    ers_den <- sapply(unique(leaves[[t]]), FUN = function(leaf)
      sapply(step, FUN = function(i)  length(y[phi.trues >= i & leaves[[t]] == leaf])))

    areas_num <- sapply(1:length(unique(leaves[[t]])), FUN = function(leaf) sum(sapply(2:length(step), FUN = function(x) s * (ers_num[x-1, leaf] + ers_num[x, leaf])/ 2 )))
    areas_den <- sapply(1:length(unique(leaves[[t]])), FUN = function(leaf) sum(sapply(2:length(step), FUN = function(x) s * (ers_den[x-1, leaf] + ers_den[x, leaf])/ 2 )))

    gammas[[t]] <- setNames(areas_num/areas_den, unique(leaves[[t]]))


    df$pseudo_res <- NULL
    strongpreds <- sapply(1:rows, FUN = function(i) strongpreds[i] + eta * gammas[[t]][names(gammas[[t]]) == leaves[[t]][names(leaves[[t]]) == i]])


    if(verbose == 1) {
      erra[t] <- sera(trues = y, preds = strongpreds, phi.trues = phi.trues)
      print(paste0("Iteration: ", t, " SERA: ", erra[t]))
    }

    t = t+1

  }

  n = nrow(test)
  m = length(stumps)

  finalpreds <- c()
  leaves_preds <- list()

  for (i in 1:m) {

    leaves_preds[[i]] <- treeClust::rpart.predict.leaves(stumps[[i]], test)
    if (is.null(names(leaves_preds[[i]]))) names(leaves_preds[[i]]) <- seq(1,n)
  }

  finalpreds <- F_0 + sapply(1:n, FUN = function(i)
    eta * sum(sapply(1:m, FUN = function(m) gammas[[m]][names(gammas[[m]]) == leaves_preds[[m]][names(leaves_preds[[m]]) == i]])))

  return(finalpreds)
}
