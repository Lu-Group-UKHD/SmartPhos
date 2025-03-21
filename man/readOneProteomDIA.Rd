% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper.R
\name{readOneProteomDIA}
\alias{readOneProteomDIA}
\title{Read and Process a Single DIA Proteomics Sample}
\usage{
readOneProteomDIA(inputTab, sampleName)
}
\arguments{
\item{inputTab}{A \code{data.table} or \code{data.frame} containing the
input data.}

\item{sampleName}{A \code{character} string specifying the sample name.}
}
\value{
A \code{data.table} containing the processed data with columns for
intensity, UniProt ID, and gene name.
}
\description{
\code{readOneProteomDIA} reads and processes data from a single DIA
proteomics sample, applying filtering and data transformation steps.
}
\details{
This function processes DIA proteomics data for a single sample by filtering
out rows with non-quantitative data, converting character values to numeric,
and renaming columns for consistency. It also ensures that each protein group
has a unique identifier.
}
