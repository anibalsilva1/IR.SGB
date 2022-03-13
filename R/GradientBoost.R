#' Determines predictions of a given dataset using the Gradient Boost for Regression using MSE as
#' an optimisation loss function.
#'
#' @param formula A formula object.
#' @param train The training dataset. An data.frame or tibble object.
#' @param test The test dataset. An data.frame or tibble object.
#' @param maxIter The maximum number of iterations.
#' @param eta Learning rate.
#' @param verbose Prints out the error across iterations (if 1).
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
#'
#' n <- nrow(NO2Emissions)
#' s <- sample(1:n, size = n*0.8)
#'
#' formula <- LNO2 ~ .
#' train <- NO2Emissions %>% slice(s)
#' test <- NO2Emissions %>% slice(-s)
#'
#' res <- GradientBoost(formula, train, test)
#' res}



GradientBoost <- function(formula,
                          train,
                          test,
                          maxIter = 200,
                          eta = 0.01,
                          verbose = 0){

  start_train_time <- Sys.time()

  rownames(train) <- 1:nrow(train)
  rownames(test) <- 1:nrow(test)

  target <- formula[[2]]

  y <- dplyr::pull(train, target)
  X <- dplyr::select(train, -target)

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


  end_train_time <- Sys.time()
  start_test_time <- Sys.time()


  while (t <= maxIter) {

    X$pseudo_res <- y - strongpreds
    X <- X %>% relocate(pseudo_res)


    resWeak <- rpart::rpart(formula = pseudo_res ~ .,
                            data = X)

    weakpreds <- predict(resWeak, X)
    stumps[[t]] <- resWeak
    gammas[t] <- sum(weakpreds * X$pseudo_res)/sum(weakpreds^2)


    strongpreds <- strongpreds + eta * gammas[t] * weakpreds

    X <- dplyr::select(X, -pseudo_res)



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

  end_test_time <- Sys.time()

  train_time <- as.numeric(difftime(end_train_time, start_train_time, units = "sec"))
  test_time <- as.numeric(difftime(end_test_time, start_test_time, units = "sec"))


  time <- c("train" = train_time, "test" = test_time)

  return(list("preds" = finalpreds,
              "time" = time))

  return(finalpreds)
}
