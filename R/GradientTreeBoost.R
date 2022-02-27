#' Determines predictions of a given dataset using the Gradient Tree Boost for Regression using
#' MSE as an optimisation loss function.
#'
#' @param formula A formula object.
#' @param train A data.frame or tibble object. The training dataset.
#' @param test A data.frame or tibble object. The test dataset.
#' @param maxIter The maximum number of iterations.
#' @param eta Learning rate.
#' @param verbose Prints out the error across iterations (if 1)
#'
#' @return A vector with predictions.
#' @export
#'
#' @examples
#' \dontrun{
#' library(IR.SGB)
#' library(dplyr)
#' n <- nrow(NO2Emissions)
#' s <- sample(1:n, size = n*0.8)
#'
#' formula <- LNO2 ~ .
#' train <- NO2Emissions %>% slice(s)
#' test <- NO2Emissions %>% slice(-s)
#'
#' preds <- GradientTreeBoost(formula, train, test)
#' preds
#' }

GradientTreeBoost <- function(formula,
                              train,
                              test,
                              maxIter = 100,
                              eta = 0.01,
                              verbose = 0){

  ### FIX ROWNAMES BECAUSE OF UNORDERED LEAF NAMES
  rownames(train) <- 1:nrow(train)
  rownames(test) <- 1:nrow(test)

  target <- formula[[2]]

  y <- dplyr::pull(train, target)
  X <- dplyr::select(train, -target)


  stumps <- list()
  gammas <- list()
  leaves <- list()

  rows = nrow(train)

  weakpreds   <- c()
  strongpreds <- c()
  pseudo_res  <- c()

  erra <- c()

  t = 1
  F_0 <- mean(y)
  strongpreds <- rep(F_0, rows)

  while (t <= maxIter) {

    X$pseudo_res <- y - strongpreds
    resWeak <- rpart::rpart(formula = pseudo_res ~ ., data = X)
    leaves[[t]] <- resWeak$where

    stumps[[t]] <- resWeak
    gammas[[t]] <- stats::setNames(sapply(unique(leaves[[t]]), FUN = function(leaf) mean(X$pseudo_res[leaves[[t]] == leaf])),
                            unique(leaves[[t]]))

    nG <- as.numeric(names(gammas[[t]]))
    nL <- as.numeric(names(leaves[[t]]))


    X$pseudo_res <- NULL

    strongpreds <- sapply(1:rows, FUN = function(i) strongpreds[i] + eta * gammas[[t]][nG == leaves[[t]][nL == i]])



    if(verbose == 1){
      erra[t] <- IRon::mse(y, strongpreds)
      print(paste0("Iteration: ", t, " mse: ", erra[t]))
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
    eta * sum(sapply(1:m,
                     FUN = function(m) gammas[[m]][names(gammas[[m]]) == leaves_preds[[m]][names(leaves_preds[[m]]) == i]])))

  return(finalpreds)
}
