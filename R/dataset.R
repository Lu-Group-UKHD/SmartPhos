# Documentation of datasets''


#' swissProt
#'
#' This is a high-quality, manually curated protein sequence database which
#' provides a high level of annotations (such as the description of the function
#' of a protein, structure of its domains, post-translational modifications,
#' variants, etc.), a minimal level of redundancy and high level of integration
#' with other databases.
#'
#' @docType data
#' @usage data(swissProt)
#' @format an object of "tbl_df" (tidy table)
#' @examples
#' data(swissProt)
#' @return A \code{data frame} or \code{tibble} containing high-level
#' annotations for manually curated proteins.
"swissProt"


#' Homo_sapien_kinase_substrate_network
#'
#' A prior knowledge database about the known kinase-phosphosite interactions
#' for Homo sapiens
#'
#' @docType data
#' @usage data(Homo_sapien_kinase_substrate_network)
#' @format a data.frame object
#' @examples
#' data(Homo_sapien_kinase_substrate_network)
#' @return A \code{data frame} containing the information about the known
#' kinase-phosphosite interactions for Homo sapiens.
"Homo_sapien_kinase_substrate_network"


#' Mus_musculus_kinase_substrate_network
#'
#' A prior knowledge database about the known kinase-phosphosite interactions
#' for Mus musculus
#'
#' @docType data
#' @usage data(Mus_musculus_kinase_substrate_network)
#' @format a data.frame object
#' @examples
#' data(Mus_musculus_kinase_substrate_network)
#' @return A \code{data frame} containing the information about the known
#' kinase-phosphosite interactions for Mus musculus.
"Mus_musculus_kinase_substrate_network"


#' dda_example
#'
#' A sample of Data Dependent Acquisition (DDA) mass spectrometry data from
#' Max Quant.
#'
#' @docType data
#' @usage data(dda_example)
#' @format a class S4 object of MultiAssayExperiment
#' @examples
#' data(dda_example)
#' @return A \code{MultiAssayExperiment} object containing a sample of DDA mass
#' spectrometry data from Max Quant.
"dda_example"


#' dia_example
#'
#' A sample of Data Independent Acquisition (DIA) mass spectrometry data from
#' Spectronaut.
#'
#' @docType data
#' @usage data(dia_example)
#' @format a class S4 object of MultiAssayExperiment
#' @examples
#' data(dia_example)
#' @return A \code{MultiAssayExperiment} object containing a sample of DIA mass
#' spectrometry data from Spectronaut.
"dia_example"
