GradientBoost <- function(formula,
                          train,
                          test,
                          weakLearner = "rpart",
                          maxIter = 200,
                          eta = 0.01,
                          verbose = 0,
                          lambda = 10,
                          maxdepth = 5){

  rownames(train) <- 1:nrow(train)
  rownames(test) <- 1:nrow(test)

  target <- formula[[2]]

  y <- pull(train, target)
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

    if(weakLearner == "lm"){

      resWeak <- do.call(weakLearner, args = list(
        formula = pseudo_res ~ .,
        data = df))
    }
    else if(weakLearner == "rpart"){

      resWeak <- do.call(weakLearner, args = list(
        formula = pseudo_res ~ .,
        data = df,
        maxdepth = maxdepth))
    }
    else if(weakLearner == "nnet"){

      resWeak <- do.call(weakLearner, args = list(
        formula = pseudo_res ~ .,
        data = df,
        size = round(sqrt(ncol(data)-1)),
        maxit = maxIter.w,
        linout =T,
        trace = F))
    }

    weakpreds <- predict(resWeak, X)
    stumps[[t]] <- resWeak
    gammas[t] <- sum(weakpreds * df$pseudo_res)/sum(weakpreds^2 + lambda)

    df$pseudo_res <- NULL
    strongpreds <- strongpreds + eta * gammas[t] * weakpreds



    if(verbose == 1){
      erra[t] <- mse(y, strongpreds)
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
