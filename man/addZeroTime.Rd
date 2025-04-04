% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/time_series_clustering.R
\name{addZeroTime}
\alias{addZeroTime}
\title{Add Zero Timepoint Data to Treatment Subset}
\usage{
addZeroTime(data, condition, treat, zeroTreat, timeRange)
}
\arguments{
\item{data}{A \code{SummarizedExperiment} object containing the experimental
data.}

\item{condition}{\code{Character} string corresponds to one of the columns
from the colData of SE object.}

\item{treat}{\code{Character} string specifying the treatment to which zero
timepoint should be added.}

\item{zeroTreat}{\code{Character} string specifying the treatment
representing the zero timepoint.}

\item{timeRange}{\code{Character} vector specifying the timepoints to include
for the treatment.}
}
\value{
A \code{SummarizedExperiment} object with the zero timepoint added to
the specified treatment's data.
}
\description{
\code{addZeroTime} adds a zero timepoint to a specific treatment's data
subset.
}
\details{
The function performs the following steps:
\enumerate{
  \item Subsets the data for the specified treatment and time range.
  \item Subsets the data for the zero timepoint of the specified zero
  treatment.
  \item Combines the assays from the treatment and zero timepoint subsets.
  \item Updates the column data to reflect the combined treatment.
  \item Returns a \code{SummarizedExperiment} object with the combined data.
}
}
\examples{
library(SummarizedExperiment)
# Load multiAssayExperiment object
data("dia_example")
# Get SummarizedExperiment object
se <- dia_example[["Phosphoproteome"]]
colData(se) <- colData(dia_example)
# Call the function
addZeroTime(se, condition = "treatment", treat = "EGF",
zeroTreat = "1stCrtl", timeRange = c("20min","40min", "6h"))

}
