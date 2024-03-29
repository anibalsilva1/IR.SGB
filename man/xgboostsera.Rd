% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{xgboostsera}
\alias{xgboostsera}
\title{Gradients Statistics of SERA for XGBoost.}
\usage{
xgboostsera(preds, dtrain)
}
\arguments{
\item{preds}{Numeric \code{vector} of predictions.}

\item{dtrain}{\code{xgb.DMatrix} object.}
}
\value{
Returns a list containing the gradient and the hessian of SERA.
}
\description{
Determines the gradient statistics of SERA to be used in
XGBoost algorithm as a custom loss function.
}
