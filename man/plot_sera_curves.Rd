% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_sera_curves.R
\name{plot_sera_curves}
\alias{plot_sera_curves}
\title{Plot SERA Curve}
\usage{
plot_sera_curves(preds, ds_name, return.err = FALSE)
}
\arguments{
\item{preds}{A \code{list} of datasets with the respective predictions.}

\item{ds_name}{A character with the dataset name.}

\item{return.err}{Returns SERA errors if TRUE. Default FALSE.}
}
\value{
SERA curves plot.
}
\description{
Plots SERA curves for several models in a dataset.
}
