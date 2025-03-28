% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataImport.R
\name{generateInputTable}
\alias{generateInputTable}
\title{Generate Input Table for Proteomic and Phosphoproteomic Analysis}
\usage{
generateInputTable(rawFolder, batchAsFolder = FALSE)
}
\arguments{
\item{rawFolder}{A \code{character} string specifying the path to the folder
containing the raw files.}

\item{batchAsFolder}{A \code{logical} value indicating whether to treat
subdirectories as separate batches. Default is \code{FALSE}.}
}
\value{
A \code{data.frame} with columns fileName, sample, searchType, batch,
and id that can be used as input for further analysis.
}
\description{
\code{generateInputTable} generates an input table for proteomic and
phosphoproteomic analysis by reading files from a specified folder.
}
\details{
The function performs the following steps:
\itemize{
  \item Optionally treats subdirectories as separate batches.
  \item Reads the summary file containing experimental information.
  \item Generates unique experimental IDs based on batch and sample names.
  \item Processes file paths for full proteome and phosphoproteome data.
  \item Creates a combined input table with file names, sample names, search
  types, batches, and IDs.
}
}
