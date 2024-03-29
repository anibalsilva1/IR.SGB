#' Gradient Tree Boosting
#'
#' @description Determines predictions of a given dataset using the Gradient Tree Boost for Regression using
#' MSE as an optimisation loss function.
#'
#' @param formula A \code{formula object}.
#' @param train A \code{data.frame} or \code{tibble} object with the training set.
#' @param test A \code{data.frame} or \code{tibble} object with the test set.
#' @param maxIter The maximum number of iterations.
#' @param eta Learning rate.
#' @param verbose Prints out the error across iterations (if 1)
#'
#' @return A numeric vector with predictions and execution time (in seconds).
#' @export
#'
#' @examples
#' \dontrun{
#'
#' library(IR.SGB)
#' library(dplyr)
#' library(rpart)
#' library(treeClust)
#'
#' n <- nrow(NO2Emissions)
#' s <- sample(1:n, size = n*0.8)
#'
#' formula <- LNO2 ~ .
#' train <- NO2Emissions %>% slice(s)
#' test <- NO2Emissions %>% slice(-s)
#'
#' res <- GradientTreeBoost(formula, train, test)
#' res
#' }

GradientTreeBoost <- function(formula,
                              train,
                              test,
                              maxIter = 100,
                              eta = 0.01,
                              verbose = 0){

  start_train_time <- Sys.time()

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

    X <- X %>% relocate(pseudo_res)

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

  end_train_time <- Sys.time()
  start_test_time <- Sys.time()

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

  end_test_time <- Sys.time()

  train_time <- as.numeric(difftime(end_train_time, start_train_time, units = "sec"))
  test_time <- as.numeric(difftime(end_test_time, start_test_time, units = "sec"))

  time <- c("train" = train_time, "test" = test_time)

  return(list("preds" = finalpreds,
              "time" = time))
}
