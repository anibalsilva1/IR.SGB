% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_diffs_bayes.R
\name{get_diffs_bayes}
\alias{get_diffs_bayes}
\title{Compares and returns the Bayes Sign Test posterior probabilities between an
oracle and the best workflows for each model.}
\usage{
get_diffs_bayes(
  res,
  bestmodels,
  oracle,
  metric = "nrmse",
  folds = 10,
  bayestest = T
)
}
\arguments{
\item{res}{A ComparisonResults object.}

\item{bestmodels}{A character vector returned from \code{\link{get_best_models_bayes}}.}

\item{oracle}{A character with the workflow name which we want to compare to \code{bestModels}.}

\item{metric}{A character with the metric name.}

\item{folds}{A value with the number of folds.}

\item{bayestest}{A boolean that specifies if we want to perform Bayes Sign Test.
If \code{False}, it only returns the average differences between the \code{oracle} and
each workflow. \code{True} by default.}
}
\value{
A list with the posterior probabilities between the \code{oracle} and
each model in \code{bestmodels} and the average differences between the \code{oracle}
and each workflow
}
\description{
Compares and returns the Bayes Sign Test posterior probabilities between an
oracle and the best workflows for each model.
}
