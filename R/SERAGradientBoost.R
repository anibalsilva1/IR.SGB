#' Gradient Boosting (SERA)
#'
#'
#' @description Determines predictions of a given dataset using the Gradient Boost for Regression using
#' SERA as an optimisation loss function.
#'
#' @param formula A \code{formula object}.
#' @param train A \code{data.frame} or \code{tibble} object with the training set.
#' @param test A \code{data.frame} or \code{tibble} object with the test set.
#' @param maxIter The maximum number of iterations.
#' @param eta Learning rate.
#' @param verbose Prints out the error across iterations (if 1).
#'
#' @return A numeric vector with predictions and execution time (in seconds).
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
#' res <- SERAGradientTreeBoost(formula, train, test)
#' res
#' }

SERAGradientBoost <- function(formula,
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
  gammas <- c()


  phi.trues <- IRon::phi(y)
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
    X <- X %>% relocate(pseudo_res)

    resWeak <- rpart::rpart(formula = pseudo_res ~ .,
                            data = X)

    weakpreds <- predict(resWeak, X)
    stumps[[t]] <- resWeak

    errs_num <- sapply(step, FUN = function(i) sum(weakpreds[phi.trues >= i]*(y[phi.trues >= i] - strongpreds[phi.trues >= i])))
    errs_den <- sapply(step, FUN = function(i) sum(weakpreds[phi.trues >= i]^2))

    areas_num <- sum(sapply(2:length(step), FUN = function(x) s * (errs_num[x-1] + errs_num[x])/ 2 ))
    areas_den <- sum(sapply(2:length(step), FUN = function(x) s * (errs_den[x-1] + errs_den[x])/ 2 ))

    gammas[t] <- areas_num/areas_den

    strongpreds <- strongpreds + eta * gammas[t] * weakpreds
    X <- dplyr::select(X, -pseudo_res)



    if(verbose == 1){
      erra[t] <- IRon::sera(trues = y, preds = strongpreds, phi.trues = phi.trues)
      print(paste0("Iteration: ", t, " SERA: ", erra[t]))
    }

    t = t+1

  }
  end_train_time <- Sys.time()
  start_test_time <- Sys.time()

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
}
