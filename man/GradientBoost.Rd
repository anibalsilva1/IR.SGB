% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GradientBoost.R
\name{GradientBoost}
\alias{GradientBoost}
\title{Determines predictions of a given dataset using the Gradient Boost for Regression using MSE as
an optimisation loss function.}
\usage{
GradientBoost(formula, train, test, maxIter = 200, eta = 0.01, verbose = 0)
}
\arguments{
\item{formula}{A formula object.}

\item{train}{The training dataset. An data.frame or tibble object.}

\item{test}{The test dataset. An data.frame or tibble object.}

\item{maxIter}{The maximum number of iterations.}

\item{eta}{Learning rate.}

\item{verbose}{Prints out the error across iterations (if 1).}
}
\value{
A numeric vector with predictions and execution time (in seconds).
}
\description{
Determines predictions of a given dataset using the Gradient Boost for Regression using MSE as
an optimisation loss function.
}
\examples{
\dontrun{

library(IR.SGB)
library(dplyr)
library(rpart)

n <- nrow(NO2Emissions)
s <- sample(1:n, size = n*0.8)

formula <- LNO2 ~ .
train <- NO2Emissions \%>\% slice(s)
test <- NO2Emissions \%>\% slice(-s)

res <- GradientBoost(formula, train, test)
res}
}
