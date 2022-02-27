#' Determines predictions of a given dataset using the Gradient Boost for Regression using MSE as
#' an optimisation loss function.
#'
#' @param formula A formula object.
#' @param train The training dataset. An data.frame or tibble object.
#' @param test The test dataset. An data.frame or tibble object.
#' @param maxIter The maximum number of iterations.
#' @param eta Learning rate.
#' @param verbose Prints out the error across iterations (if 1).
#' @param lambda Ridge regression parameter.
#' @param maxdepth Max depth of a tree.
#'
#' @return A numeric vector with predictions.
#' @export
#'
#' @examples
#' \dontrun{
#' library(IR.SGB)
#' library(dplyr)
#'
#' n <- nrow(NO2Emissions)
#' s <- sample(1:n, size = n*0.8)
#'
#' formula <- LNO2 ~ .
#' train <- NO2Emissions %>% slice(s)
#' test <- NO2Emissions %>% slice(-s)
#'
#' preds <- GradientBoost(formula, train, test)
#' preds}



GradientBoost <- function(formula,
                          train,
                          test,
                          maxIter = 200,
                          eta = 0.01,
                          verbose = 0,
                          lambda = 10,
                          maxdepth = 5){

  rownames(train) <- 1:nrow(train)
  rownames(test) <- 1:nrow(test)

  target <- formula[[2]]

  y <- dplyr::pull(train, target)
  X <- dplyr::select(train, -target)
  df <- X

  stumps <- list()
  gammas <- c()

  rows = nrow(train)

  weakpreds   <- c()
  strongpreds <- c()
  pseudo_res  <- c()
  erra <- c()

  t = 1
  F_0 <- mean(y)
  strongpreds <- rep(F_0, rows)


  while (t <= maxIter) {

    df$pseudo_res <- y - strongpreds


    resWeak <- rpart::rpart(formula = pseudo_res ~ .,
                            data = df,
                            maxdepth = maxdepth)

    weakpreds <- predict(resWeak, X)
    stumps[[t]] <- resWeak
    gammas[t] <- sum(weakpreds * df$pseudo_res)/sum(weakpreds^2 + lambda)

    df$pseudo_res <- NULL
    strongpreds <- strongpreds + eta * gammas[t] * weakpreds



    if(verbose == 1){
      erra[t] <- IRon::mse(y, strongpreds)
      print(paste0("Iteration: ", t, " mse: ", erra[t]))
    }
    t = t+1

  }

  m = nrow(test)
  n = length(stumps)
  preds <- matrix(0, nrow = m, ncol = n)
  finalpreds <- c()

  for (i in 1:n) {

    preds[, i] <- predict(stumps[[i]], test)

  }
  finalpreds <- F_0 + sapply(1:m, FUN = function(i) eta * sum(gammas * preds[i ,]))

  return(finalpreds)
}
