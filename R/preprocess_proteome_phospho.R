#' @name getOneSymbol
#'
#' @title Extract the Last Gene Symbol from a Semicolon-Separated List
#'
#' @description
#' \code{getOneSymbol} extracts the last gene symbol from a semicolon-separated list of gene symbols.
#'
#' @param Gene A \code{character} vector where each element is a semicolon-separated list of gene symbols.
#'
#' @return A \code{character} vector containing the last gene symbol from each element of the input vector.
#'
#' @details
#' This function processes a character vector where each element consists of gene symbols separated by semicolons. It splits each element by semicolons and extracts the last gene symbol from the resulting list. The output is a character vector of these last gene symbols.
#'
#' @importFrom stringr str_split
getOneSymbol <- function(Gene) {
  # Apply a function to each element in the Gene vector
  outStr <- sapply(Gene, function(x) {
    # Split the string by semicolons
    sp <- str_split(x, ";")[[1]]
    # Return the last element in the split string
    sp[length(sp)]
  })
  # Remove names from the resulting vector
  names(outStr) <- NULL
  return(outStr)
}


#' @name preprocessProteome
#'
#' @title Preprocess Proteome Data
#'
#' @description
#' \code{preprocessProteome} preprocesses proteome data stored in a \code{SummarizedExperiment} object by performing filtering, transformation, normalization, imputation, and batch effect removal.
#'
#' @param seData A \code{SummarizedExperiment} object containing proteome data.
#' @param filterList A \code{list} of filters to apply on the samples. Default is \code{NULL}.
#' @param missCut \code{Numeric} value specifying the missing value cutoff percentage for filtering features. Default is 50.
#' @param transform \code{Character} string specifying the transformation method ("log2", "vst", "none"). Default is "log2".
#' @param normalize \code{Logical} value indicating whether to normalize the data. Default is \code{FALSE}.
#' @param getPP \code{Logical} value indicating whether to retrieve PP samples. Default is \code{FALSE}.
#' @param removeOutlier \code{Character} vector of samples to be removed as outliers. Default is \code{NULL}.
#' @param impute \code{Character} string specifying the imputation method ("QRILC", "MLE", "bpca", "missForest", "MinDet", "none"). Default is "none".
#' @param batch \code{Character} vector specifying batch effects to remove. Default is \code{NULL}.
#' @param verbose \code{Logical} value indicating whether to print detailed information. Default is \code{FALSE}.
#' @param scaleFactorTab \code{Data frame} containing scale factors for normalization. Default is \code{NULL}.
#'
#' @return A \code{SummarizedExperiment} object with preprocessed proteome data.
#'
#' @examples
#' library(SummarizedExperiment)
#' # Load multiAssayExperiment object
#' data("dia_example")
#' # Get SummarizedExperiment object
#' se <- dia_example[["Proteome"]]
#' colData(se) <- colData(dia_example)
#' # Call the function
#' preprocessProteome(seData = se, normalize = TRUE, impute = "QRILC")
#'
#' @importFrom SummarizedExperiment colData rowData assay assays
#' @importFrom dplyr filter mutate
#' @importFrom tidyr pivot_longer
#' @importFrom vsn justvsn
#' @importFrom DEP impute
#' @importFrom limma removeBatchEffect
#' @importFrom missForest missForest
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG registerDoRNG
#'
#' @export
preprocessProteome <- function(seData, filterList = NULL, missCut = 50,
                               transform = "log2", normalize = FALSE, getPP = FALSE,
                               removeOutlier = NULL, impute = "none", batch = NULL,
                               verbose = FALSE, scaleFactorTab = NULL) {

  # Retrieve desired sample type
  if (getPP) {
    # Retrieve PP samples if specified
    fpe <- seData[,seData$sampleType %in% c("Phospho", "PP")]
    colData(fpe) <- colData(seData)[colnames(fpe),]
  } else {
    # Otherwise, retrieve FullProteome samples
    fpe <- seData[,seData$sampleType %in% c("FullProteome", "FP")]
    colData(fpe) <- colData(seData)[colnames(fpe),]
  }

  # Remove specified outliers
  if (length(removeOutlier) > 0) {
    if (length(removeOutlier) > 1) {
      for (i in removeOutlier) {
        fpe <- fpe[, !grepl(i, fpe$sample)]
      }
    }
    else {
      fpe <- fpe[, !grepl(removeOutlier, fpe$sample)]
    }
  }

  # Apply specified filters
  if (!is.null(filterList)) {
    for (n in names(filterList)) {
      fpe <- fpe[,fpe[[n]] %in% filterList[[n]]]
    }
  }

  # Rename columns to sample names
  colnames(fpe) <- fpe$sample

  # Process gene names and remove features without symbols
  rowData(fpe)$Gene <- getOneSymbol(rowData(fpe)$Gene)
  fpe <- fpe[!rowData(fpe)$Gene %in% c(NA,""),]

  # Filter features based on missing values
  countMat <- assay(fpe)
  missPer <- rowSums(is.na(countMat))/ncol(countMat)*100
  fpeSub <- fpe[missPer < missCut,]

  # Apply transformation and normalization
  if (transform=="log2") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(fpeSub) <- medianNorm(log2(assay(fpeSub)))
      } else {
        assay(fpeSub) <- log2(t(t(assay(fpeSub))/scaleFactorTab[match(paste0(fpeSub$sample),scaleFactorTab$sample),]$scaleFactor))
      }
    } else {
      assay(fpeSub) <- log2(assay(fpeSub))
    }
  } else if (transform == "vst") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(fpeSub) <- vsn::justvsn(assay(fpeSub))
      } else {
        normMat <- t(t(assay(fpeSub))/scaleFactorTab[match(paste0(fpeSub$sample),scaleFactorTab$sample),]$scaleFactor)
        assay(fpeSub) <- vsn::justvsn(normMat, calib="none")
      }
    } else {
      assay(fpeSub) <- vsn::justvsn(assay(fpeSub), calib="none")
    }
  } else if (transform == "none") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(fpeSub) <- medianNorm(assay(fpeSub))
      } else {
        assay(fpeSub) <- t(t(assay(fpeSub))/scaleFactorTab[match(paste0(fpeSub$sample),scaleFactorTab$sample),]$scaleFactor)
      }
    } else {
      assay(fpeSub) <- assay(fpeSub)
    }
  }

  # Impute missing values
  if (impute != "none") {
    rowData(fpeSub)$name <- rowData(fpeSub)$UniprotID
    rowData(fpeSub)$ID <- rowData(fpeSub)$UniprotID
    if (impute == "missForest") {
        doParallel::registerDoParallel(cores = 6)  # set based on number of CPU cores
        doRNG::registerDoRNG(seed = 123)
        mf <- missForest::missForest(t(assay(fpeSub)), parallelize = "forests", maxiter = 2, ntree = 50)
        imp <- t(mf$ximp)
    }
    else {
      imp <- DEP::impute(fpeSub, fun = impute)
    }
    assays(fpeSub)[["imputed"]] <- assay(imp)
    rowData(fpeSub)$name <- NULL
    rowData(fpeSub)$ID <- NULL
    if (verbose) {
      # Show number of samples and features
      message("Number of proteins and samples:")
      print(dim(fpeSub))
    }
  }

  # Remove batch effects if specified
  if(!is.null(batch)) {
    if(length(batch) == 1) {
      remBatchImp <- limma::removeBatchEffect(assays(fpeSub)[["imputed"]],
                                              batch = colData(fpeSub)[,batch])
      remBatch <- limma::removeBatchEffect(assay(fpeSub),
                                           batch = colData(fpeSub)[,batch])
    }
    else {
      remBatchImp <- limma::removeBatchEffect(assays(fpeSub)[["imputed"]],
                                              batch = colData(fpeSub)[,batch[1]],
                                              batch2 = colData(fpeSub)[,batch[2]])
      remBatch <- limma::removeBatchEffect(assay(fpeSub),
                                           batch = colData(fpeSub)[,batch[1]],
                                           batch2 = colData(fpeSub)[,batch[2]])
    }
    assays(fpeSub)[["imputed"]] <- assay(remBatchImp)
    assay(fpeSub) <- assay(remBatch)
  }

  return(fpeSub)
}

#' @name preprocessPhos
#'
#' @title Preprocess Phosphoproteome Data
#'
#' @description
#' \code{preprocessPhos} preprocesses phosphoproteome data stored in a \code{SummarizedExperiment} object by performing filtering, transformation, normalization, imputation, and batch effect removal.
#'
#' @param seData A \code{SummarizedExperiment} object containing phosphoproteome data.
#' @param filterList A \code{list} of filters to apply on the samples. Default is \code{NULL}.
#' @param missCut \code{Numeric} value specifying the missing value cutoff percentage for filtering features. Default is 50.
#' @param transform \code{Character} string specifying the transformation method ("log2", "vst", "none"). Default is "log2".
#' @param normalize \code{Logical} value indicating whether to normalize the data. Default is \code{FALSE}.
#' @param getFP \code{Logical} value indicating whether to retrieve FP samples. Default is \code{FALSE}.
#' @param removeOutlier \code{Character} vector of samples to be removed as outliers. Default is \code{NULL}.
#' @param assayName \code{Character} string specifying the assay name in the SummarizedExperiment object. Default is \code{NULL}.
#' @param batch \code{Character} vector specifying batch effects to remove. Default is \code{NULL}.
#' @param scaleFactorTab \code{Data frame} containing scale factors for normalization. Default is \code{NULL}.
#' @param impute \code{Character} string specifying the imputation method ("QRILC", "MLE", "bpca", "missForest", "MinDet", "none"). Default is "none".
#' @param verbose \code{Logical} value indicating whether to print detailed information. Default is \code{FALSE}.
#'
#' @return A \code{SummarizedExperiment} object with preprocessed phosphoproteome data.
#'
#' @examples
#' library(SummarizedExperiment)
#' # Load multiAssayExperiment object
#' data("dia_example")
#' # Get SummarizedExperiment object
#' se <- dia_example[["Phosphoproteome"]]
#' colData(se) <- colData(dia_example)
#' # Call the function
#' preprocessPhos(seData = se, normalize = TRUE, impute = "QRILC")
#'
#' @importFrom SummarizedExperiment colData rowData assay assays
#' @importFrom dplyr filter mutate
#' @importFrom tidyr pivot_longer
#' @importFrom DEP impute
#' @importFrom vsn justvsn
#' @importFrom limma removeBatchEffect
#' @importFrom missForest missForest
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG registerDoRNG
#'
#' @export
preprocessPhos <- function(seData, filterList = NULL, missCut = 50,
                           transform="log2", normalize = FALSE, getFP = FALSE,
                           removeOutlier = NULL, assayName = NULL, batch = NULL,
                           scaleFactorTab = NULL, impute = "none", verbose = FALSE) {

  # Retrieve the desired sample type or specified assay
  if (is.null(assayName)) {
    if (getFP) {
      # Retrieve FullProteome samples if specified
      ppe <- seData[,seData$sampleType %in% c("FullProteome", "FP")]
      colData(ppe) <- colData(seData)[colnames(ppe),]
    } else {
      # Otherwise, retrieve Phospho samples
      ppe <- seData[,seData$sampleType %in% c("Phospho", "PP")]
      colData(ppe) <- colData(seData)[colnames(ppe),]
    }
  } else {
    ppe <- seData[[assayName]]
    colData(ppe) <- colData(seData[,colnames(ppe)])
  }

  # Remove specified outliers
  if (length(removeOutlier) > 0) {
    if (length(removeOutlier) > 1) {
      for (i in removeOutlier) {
        ppe <- ppe[, !grepl(i, ppe$sample)]
      }
    }
    else {
      ppe <- ppe[, !grepl(removeOutlier, ppe$sample)]
    }
  }

  # Apply specified filters
  if (!is.null(filterList)) {
    for (n in names(filterList)) {
      ppe <- ppe[,ppe[[n]] %in% filterList[[n]]]
    }
  }

  # Rename columns to sample names
  colnames(ppe) <- ppe$sample
  # Get last gene name
  rowData(ppe)$Gene <- getOneSymbol(rowData(ppe)$Gene)
  # Get last phosphorylation site
  rowData(ppe)$Residue <- getOneSymbol(rowData(ppe)$Residue)
  rowData(ppe)$Position <- getOneSymbol(rowData(ppe)$Position)
  # Remove features without gene symbols
  ppe <- ppe[!rowData(ppe)$Gene %in% c(NA,""),]
  # Rename phosphorylation sites
  rowData(ppe)$site <- paste0(rowData(ppe)$Gene,"_",rowData(ppe)$Residue,rowData(ppe)$Position)
  # Filter features based on missing values
  countMat <- assay(ppe)
  missPer <- rowSums(is.na(countMat))/ncol(countMat)*100
  ppeSub <- ppe[missPer < missCut,]

  # Apply transformation and normalization
  if (transform=="log2") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(ppeSub) <- medianNorm(log2(assay(ppeSub)))
      } else {
        assay(ppeSub) <- log2(t(t(assay(ppeSub))/scaleFactorTab[match(paste0(ppeSub$sample),scaleFactorTab$sample),]$scaleFactor))
      }
    } else {
      assay(ppeSub) <- log2(assay(ppeSub))
    }
  } else if (transform == "vst") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(ppeSub) <- vsn::justvsn(assay(ppeSub))
      } else {
        normMat <- t(t(assay(ppeSub))/scaleFactorTab[match(paste0(ppeSub$sample),scaleFactorTab$sample),]$scaleFactor)
        assay(ppeSub) <- vsn::justvsn(normMat, calib="none")
      }
    } else {
      assay(ppeSub) <- vsn::justvsn(assay(ppeSub), calib="none")
    }
  } else if (transform == "none") {
    if (normalize) {
      if (is.null(scaleFactorTab)) {
        assay(ppeSub) <- medianNorm(assay(ppeSub))
      } else {
        assay(ppeSub) <- t(t(assay(ppeSub))/scaleFactorTab[match(paste0(ppeSub$sample),scaleFactorTab$sample),]$scaleFactor)
      }
    } else {
      assay(ppeSub) <- assay(ppeSub)
    }
  }

  # Impute missing values
  if (impute != "none") {
    rowData(ppeSub)$name <- rowData(ppeSub)$site
    rowData(ppeSub)$ID <- rowData(ppeSub)$site
    if (impute == "missForest") {
        doParallel::registerDoParallel(cores = 6)  # Set number of CPU cores
        doRNG::registerDoRNG(seed = 123)
        mf <- missForest::missForest(t(assay(ppeSub)), parallelize = "forests", maxiter = 2, ntree = 50)
        imp <- t(mf$ximp)
    }
    else {
      imp <- DEP::impute(ppeSub, fun = impute)
    }
    assays(ppeSub)[["imputed"]] <- assay(imp)
    rowData(ppeSub)$name <- NULL
    rowData(ppeSub)$ID <- NULL
    # Show number of samples and features
    if (verbose) {
      message("Number of proteins and samples:")
      print(dim(ppeSub))
    }
  }

  # Remove batch effects if specified
  if(!is.null(batch)) {
    if(length(batch) == 1) {
      remBatchImp <- limma::removeBatchEffect(assays(ppeSub)[["imputed"]],
                                              batch = colData(ppeSub)[,batch])
      remBatch <- limma::removeBatchEffect(assay(ppeSub),
                                           batch = colData(ppeSub)[,batch])
    }
    else {
      remBatchImp <- limma::removeBatchEffect(assays(ppeSub)[["imputed"]],
                                              batch = colData(ppeSub)[,batch[1]],
                                              batch2 = colData(ppeSub)[,batch[2]])
      remBatch <- limma::removeBatchEffect(assay(ppeSub),
                                           batch = colData(ppeSub)[,batch[1]],
                                           batch2 = colData(ppeSub)[,batch[2]])
    }
    assays(ppeSub)[["imputed"]] <- assay(remBatchImp)
    assay(ppeSub) <- assay(remBatch)
  }

  return(ppeSub)
}
