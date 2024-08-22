#' @name medianNorm
#' 
#' @title Normalize a Matrix Using Median or Mean
#'
#' @description
#' `medianNorm` normalizes the columns of a matrix by either the median or the mean.
#'
#' @param x A numeric matrix to be normalized.
#' @param method A character string specifying the normalization method. Options are `"median"` or `"mean"`. Default is `"median"`.
#'
#' @return A numeric matrix with normalized columns.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item If the `method` is `"median"`, it calculates the median of each column and adjusts by the overall median of these medians.
#'   \item If the `method` is `"mean"`, it calculates the mean of each column and adjusts by the overall mean of these means.
#'   \item It constructs a matrix of these adjusted values and subtracts it from the original matrix to normalize the columns.
#' }
#'
#' @examples
#' # Example usage:
#' x <- matrix(rnorm(20), nrow=5, ncol=4)
#' normalized_x <- medianNorm(x, method = "median")
#' print(normalized_x)
#'
#' @importFrom matrixStats colMedians
#' @export
medianNorm <- function(x, method = "median") {
  if (method == "median") {
    # Calculate the median of each column, ignoring NA values
    mVal <- matrixStats::colMedians(x, na.rm = TRUE)
    # Adjust by the overall median of these medians
    mVal <- mVal - median(mVal, na.rm = TRUE)
  } else if (method == "mean") {
    # Calculate the mean of each column, ignoring NA values
    mVal <- colMeans(x, na.rm = TRUE)
    # Adjust by the overall mean of these means
    mVal <- mVal - mean(mVal, na.rm = TRUE)
  }
  mMat <- matrix(rep(mVal, each = nrow(x)), ncol =ncol(x))
  return(x-mMat)
}

# function for performing normalization of FP and PP samples

#' @name performCombinedNormalization
#' 
#' @title Perform Combined Normalization on MultiAssayExperiment Data
#'
#' @description
#' `performCombinedNormalization` performs combined normalization on proteome and phosphoproteome data from a MultiAssayExperiment object.
#'
#' @param maeData A MultiAssayExperiment object containing proteome and phosphoproteome data.
#'
#' @return A numeric matrix with normalized and log2-transformed data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Extracts the count matrices for Full Proteome (FP) samples.
#'   \item Combines the proteome and phosphoproteome data into a single matrix.
#'   \item Removes rows with all NA values.
#'   \item Performs median normalization and log2 transformation on the combined matrix.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming maeData is a MultiAssayExperiment object with appropriate data
#' # normalized_data <- performCombinedNormalization(maeData)
#' # print(normalized_data)
#'
#' @importFrom MultiAssayExperiment assay
#' @export
performCombinedNormalization <- function(maeData) {
  
  # get count matrix from FP (Full Proteome) samples
  setFP <- maeData[,maeData$sampleType %in% c("FullProteome", "FP")]
  protFP <- assay(setFP[["Proteome"]])
  phosFP <- assay(setFP[["Phosphoproteome"]])
  # Combine proteome and phosphoproteome data into a single matrix
  comFP <- rbind(protFP, phosFP)
  comFP <- comFP[rowSums(!is.na(comFP))>0,]
  
  # perform median normalization and log2 transformation
  comFP.norm <- medianNorm(log2(comFP))
  
  return(comFP.norm)
}


#' @name getRatioMatrix
#' 
#' @title Get Ratio Matrix of Phosphoproteome Data
#'
#' @description
#' `getRatioMatrix` calculates the ratio matrix of phosphoproteome data from a MultiAssayExperiment object.
#'
#' @param maeData A MultiAssayExperiment object containing phosphoproteome and full proteome data.
#' @param normalization A logical value indicating whether to perform normalization. Default is `FALSE`.
#' @param getAdjustedPP A logical value indicating whether to use adjusted phosphoproteome data. Default is `FALSE`.
#'
#' @return A numeric matrix representing the ratio of intensity of PP (phosphoproteome) data to FP (full proteome) data.
#'
#' @examples
#' # Example usage:
#' # Assuming maeData is a MultiAssayExperiment object with appropriate data
#' # ratio_matrix <- getRatioMatrix(maeData, normalization = TRUE, getAdjustedPP = TRUE)
#' # print(ratio_matrix)
#'
#' @importFrom MultiAssayExperiment assay assays
#' @export
getRatioMatrix <- function(maeData, normalization = FALSE, getAdjustedPP = FALSE) {
  
  # Ensure the normalization parameter is logical
  stopifnot(is.logical(normalization))
  
  # Extract and log-transform phosphoproteome data based on getAdjustedPP flag
  if (!getAdjustedPP) {
    phosPP <- log2(assay(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]]))
  } else {
    phosPP <- log2(assays(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]])[["Intensity_adjusted"]])
  }
  
  # Extract and log-transform full proteome data, with optional normalization
  if (!normalization) {
    phosFP <- log2(assay(maeData[,maeData$sampleType %in% c("FullProteome", "FP")][["Phosphoproteome"]]))
  } 
  else  {
    phosFP <- performCombinedNormalization(maeData)
  }
  
  # Use the sample name without prefix as column name
  colnames(phosFP) <- maeData[,maeData$sampleType %in% c("Phospho", "PP")]$sampleName
  colnames(phosPP) <- maeData[,maeData$sampleType %in% c("FullProteome", "FP")]$sampleName
  
  # Calculate the ratio matrix
  allSmp <- intersect(colnames(phosFP), colnames(phosPP))
  allRow <- intersect(rownames(phosFP), rownames(phosPP))
  ratioMat <- phosPP[allRow, allSmp] - phosFP[allRow, allSmp]
  ratioMat <- ratioMat[rowSums(!is.na(ratioMat)) >0,]
  
  return(ratioMat)
  
}


# plot the log ration of PP/FP intensities

#' @name plotLogRatio
#' 
#' @title Plot Log Ratio of PP/FP (Phosphoproteome to Full Proteome) intensities
#'
#' @description
#' `plotLogRatio` generates a boxplot of the log2 ratio of intensities of phosphoproteome to full proteome data from a MultiAssayExperiment object.
#'
#' @param maeData A MultiAssayExperiment object containing phosphoproteome and full proteome data.
#' @param normalization A logical value indicating whether to perform normalization. Default is `FALSE`.
#'
#' @return A ggplot object representing the boxplot of the log2 ratios.
#'
#' @examples
#' # Example usage:
#' # Assuming maeData is a MultiAssayExperiment object with appropriate data
#' # plot <- plotLogRatio(maeData, normalization = TRUE)
#' # print(plot)
#'
#' @importFrom MultiAssayExperiment assay
#' @importFrom matrixStats colMedians
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr filter mutate
#' @importFrom ggplot2 ggplot aes geom_boxplot ggtitle xlab ylab geom_hline theme element_text
#' @export 
plotLogRatio <- function(maeData, normalization = FALSE) {
  # Calculate the ratio matrix of phosphoproteome to full proteome data
  ratioMat <- getRatioMatrix(maeData, normalization)
  # Extract and log-transform phosphoproteome data
  phosPP <- log2(assay(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]]))
  medianPP <- colMedians(phosPP,na.rm = TRUE)
  names(medianPP) <- maeData[,maeData$sampleType %in% c("Phospho", "PP")]$sampleName
  
  # Create a table for plotting
  plotTab <- as_tibble(ratioMat, rownames = "feature") %>%
    pivot_longer(-feature) %>%
    filter(!is.na(value)) %>%
    mutate(medianPP = medianPP[name])
  
  # Generate a ggplot boxplot of the log2 ratios
  ggplot(plotTab, aes(x=name, y=value)) +
    geom_boxplot(aes(fill = medianPP)) +
    ggtitle("Boxplot of Phospho/FullProteome Ratio") +
    xlab("sample") +
    ylab("log2(ratio)") +
    geom_hline(yintercept = median(median(plotTab$value, na.rm=TRUE)), linetype = "dashed", color = "red") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
}

# check the PP/FP ratio matrix and remove feature that do not meet requirements

#' @name checkRatioMat
#' 
#' @title Check the PP/FP ratio matrix and remove feature that do not meet requirements
#'
#' @description
#' `checkRatioMat` checks the ratio matrix for samples that do not have sufficient overlap of phospho-peptides between enriched (PP) and unenriched (FP) samples.
#'
#' @param ratioMat A numeric matrix representing the ratio of phosphoproteome data to full proteome data.
#' @param minOverlap An integer specifying the minimum number of overlapping peptides required between samples. Default is `3`.
#'
#' @return A character vector of sample names that do not meet the overlap criteria.
#'
#' @examples
#' # Example usage:
#' # Assuming ratioMat is a matrix with appropriate data
#' # excluded_samples <- checkRatioMat(ratioMat, minOverlap = 3)
#' # print(excluded_samples)
#'
#' @export
checkRatioMat <- function(ratioMat, minOverlap = 3) {
  # Initialize a list to keep track of excluded samples
  excludeSampleList <- c()
  
  # Identify samples that don't have any phospho sites detect in both FP and PP samples
  noOverSmp <- colnames(ratioMat)[colSums(!is.na(ratioMat))==0]
  if (length(noOverSmp) >0) {
    warning(paste0("Below samples don't have phopho-peptides detected in both enriched (PP) and unenriched (FP) samples and therefore adjusting factor will set to 0 (no adjustment) for them:\n",
                   paste0(noOverSmp, collapse = ", ")))
    excludeSampleList <- c(excludeSampleList, noOverSmp)
  }
  
  # Remove the identified samples from the ratio matrix
  ratioMat <- ratioMat[, !colnames(ratioMat) %in% noOverSmp]
  
  # Check for samples that do not have enough peptide overlap with other samples
  pairOverlap <- sapply(colnames(ratioMat), function(n) {
    subMat <- ratioMat[!is.na(ratioMat[,n]),]
    minOver <- min(colSums(!is.na(subMat)))
  })
  
  tooFewOverlap <- colnames(ratioMat)[pairOverlap < minOverlap]
  
  if (length(tooFewOverlap) >0) {
    warning(paste0("Below samples don't enough number of overlapped phopho-peptides with other samples and therefore adjusting factor will set to 0 (no adjustment) for them:\n",
                   paste0(tooFewOverlap, collapse = ", ")))
    excludeSampleList <- c(excludeSampleList, tooFewOverlap)
  }
  
  return(excludeSampleList)
}



#' @name runPhosphoAdjustment
#' 
#' @title Run Phospho Adjustment
#'
#' @description
#' `runPhosphoAdjustment` performs phospho adjustment on a MultiAssayExperiment object to normalize the phosphoproteome data.
#'
#' @param maeData A MultiAssayExperiment object containing phosphoproteome and full proteome data.
#' @param normalization A logical value indicating whether to perform normalization. Default is `FALSE`.
#' @param minOverlap An integer specifying the minimum number of overlapping peptides required between samples. Default is `3`.
#' @param completeness A numeric value indicating the required completeness of data for features to be included. Default is `0`.
#' @param ncore An integer specifying the number of cores to use for parallel processing. Default is `1`.
#'
#' @return A MultiAssayExperiment object with adjusted phosphoproteome data.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Defines an optimization function to minimize the sum of squared differences between pairs of samples.
#'   \item Calculates the ratio matrix of phosphoproteome to full proteome data.
#'   \item Subsets features based on completeness criteria.
#'   \item Performs a sanity check to identify and exclude problematic samples.
#'   \item Sets initial values for the adjustment factor based on column medians.
#'   \item Estimates the adjustment factor using parallel optimization.
#'   \item Adjusts the phosphoproteome measurements using the estimated adjustment factor.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming maeData is a MultiAssayExperiment object with appropriate data
#' # adjusted_maeData <- runPhosphoAdjustment(maeData, normalization = TRUE, minOverlap = 3, completeness = 0.8, ncore = 2)
#' # print(adjusted_maeData)
#'
#' @importFrom MultiAssayExperiment assay assays colData
#' @importFrom matrixStats colMedians
#' @importFrom stats optim
#' @importFrom parallel makeCluster setDefaultCluster stopCluster
#' @importFrom utils combn
#' @importFrom optimParallel optimParallel
#' @export
runPhosphoAdjustment <- function(maeData, normalization = FALSE, minOverlap = 3, completeness = 0, ncore = 1 ) {
  
  # Function to opitmize
  esFun <- function(par, data) {
    comPair <- utils::combn(seq(length(par)), 2)
    sum(((data[comPair[1, ],] + par[comPair[1, ]]) - (data[comPair[2, ],] + par[comPair[2, ]]))^2/rowSums(!is.na(data[comPair[1,],] + data[comPair[2,],])), na.rm = TRUE)
  }
  
  # Get PP/FP ratio matrix
  ratioMat <- getRatioMatrix(maeData, normalization = normalization)
  adjFac <- structure(rep(0, length.out = ncol(ratioMat)), names = colnames(ratioMat))
  
  # Subset features according to completeness in the ratio matrix
  ratioMat <- ratioMat[rowSums(!is.na(ratioMat))/ncol(ratioMat) >= completeness,]
  
  # Sanity check to see if any samples need to be excluded
  excList <- checkRatioMat(ratioMat)
  ratioMat <- ratioMat[, !colnames(ratioMat) %in% excList]
  
  # Set an initial value for B based on col medians of ratioMat, may increase search speed
  colMed <- apply(ratioMat,2, median, na.rm = TRUE)
  iniPar <- median(colMed) - colMed
  
  # Estimating adjusting factor
  cl <- makeCluster(ncore)
  setDefaultCluster(cl = cl)
  optRes <- optimParallel(par=iniPar, fn=esFun, data=t(ratioMat))
  stopCluster(cl)
  
  # Add adjusting factor to sample annotation  
  adjFac[names(optRes$par)] <- optRes$par
  ppName <- colnames(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]])
  adjFac <- structure(adjFac[maeData[,ppName]$sampleName], names = ppName)
  maeData$adjustFactorPP <- unname(adjFac[match(rownames(colData(maeData)),names(adjFac))])
  # Adjust phospho measurement on PP samples
  phosMat <- assay(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]])
  phosMat <- t(t(phosMat)*(2^adjFac))
  assays(maeData[["Phosphoproteome"]])[["Intensity_adjusted"]] <- assays(maeData[["Phosphoproteome"]])[["Intensity"]]
  assays(maeData[["Phosphoproteome"]])[["Intensity_adjusted"]][,colnames(phosMat)] <- phosMat
  
  return(maeData)
}


#' @name plotAdjustmentResults
#' 
#' @title Plot Adjustment Results
#'
#' @description
#' `plotAdjustmentResults` generates plots to visualize the results of phosphoproteome adjustment.
#'
#' @param maeData A MultiAssayExperiment object containing phosphoproteome and full proteome data.
#' @param normalization A logical value indicating whether normalization was performed. Default is `FALSE`.
#'
#' @return A list containing:
#' \item{ratioTrendPlot}{A ggplot object showing the line plot of PP/FP ratios for features present in all samples.}
#' \item{ratioBoxplot}{A ggplot object showing the box plot of PP/FP ratios before and after adjustment.}
#' \item{ppBoxplot}{A ggplot object showing the box plot of phosphorylation intensities in PP samples before and after adjustment.}
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if the adjustment factor is present in the sample annotation.
#'   \item Calculates the ratio matrix before and after adjustment.
#'   \item Creates a trend line plot for features present in all samples.
#'   \item Creates box plots of the PP/FP ratios and phosphorylation intensities before and after adjustment.
#' }
#'
#' @examples
#' # Example usage:
#' # Assuming maeData is a MultiAssayExperiment object with appropriate data
#' # plots <- plotAdjustmentResults(maeData, normalization = TRUE)
#' # print(plots$ratioTrendPlot)
#' # print(plots$ratioBoxplot)
#' # print(plots$ppBoxplot)
#'
#' @importFrom MultiAssayExperiment assay assays colData
#' @importFrom dplyr bind_rows filter group_by mutate summarise
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 aes facet_wrap geom_boxplot geom_hline geom_line geom_point ggplot ggtitle theme xlab ylab
#' @export
plotAdjustmentResults <- function(maeData, normalization = FALSE) {
  # Check if the adjustment factor has been applied
  if (!"adjustFactorPP" %in% colnames(colData(maeData))) {
    stop("Phosphorylation measurments have not been adjusted yet. Please perform normalization adjustment using calcAdjustFacotr function first")
  }
  
  # Visualize precursors before and after adjustment
  ratioMat.ori <- getRatioMatrix(maeData, normalization = normalization, getAdjustedPP = FALSE)
  ratioMat.adj <- getRatioMatrix(maeData, normalization = normalization, getAdjustedPP = TRUE)
  ratioPlotTab <- bind_rows(pivot_longer(as_tibble(ratioMat.ori, rownames = "id"), -id, names_to = "sample", values_to = "ratio") %>% mutate(adjustment = "before adjustment"),
                            pivot_longer(as_tibble(ratioMat.adj, rownames = "id"), -id, names_to = "sample", values_to = "ratio") %>% mutate(adjustment = "after adjustment")) %>%
    mutate(adjustment = factor(adjustment, levels = c("before adjustment","after adjustment"))) %>%
    filter(!is.na(ratio))
  
  # For precursors present in all samples
  featureComplete <- rownames(ratioMat.ori)[complete.cases(ratioMat.ori)]
  if (!length(featureComplete) >0) {
    warning("No feature (PP/FP ratio) has been detected in all samples. Ratio trend line will not be generated")
    ratioTrendPlot <- NULL 
  } else {
    ratioPlotTab.complete <- filter(ratioPlotTab, id %in% featureComplete)
    ratioTrendPlot <- ggplot(ratioPlotTab.complete, aes(x=sample, y=ratio, color = adjustment)) +
      geom_line(aes(group =id), linetype = "dashed") +
      geom_point() +
      facet_wrap(~adjustment, ncol=1) +
      ggtitle("Line plot of PP/FP ratios available for all samples") + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            legend.position = "none") +
      xlab("") + ylab("log2(PP/FP) ratio") 
  }
  
  # For ratio box plots
  medTab <- group_by(ratioPlotTab, adjustment) %>%
    summarise(medVal = median(ratio, na.rm=TRUE))
  ratioBoxplot <- ggplot(ratioPlotTab, aes(x=sample, y=ratio, fill = adjustment)) +
    geom_boxplot() +
    geom_hline(data=medTab, aes(yintercept = medVal), linetype = "dashed", color = "red") +
    facet_wrap(~adjustment, ncol=1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none") +
    xlab("") + ylab("log2(PP/FP) ratio") +
    ggtitle("Box plot of all PP/FP ratios")  
  
  
  # For phosphorylation measurements of PP samples before and after adjustment
  ppMat.adj <- assays(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]])[["Intensity_adjusted"]] 
  ppMat.ori <- assays(maeData[,maeData$sampleType %in% c("Phospho", "PP")][["Phosphoproteome"]])[["Intensity"]] 
  ppPlotTab <- bind_rows(pivot_longer(as_tibble(ppMat.ori, rownames = "id"), -id, names_to = "sample", values_to = "count") %>% mutate(adjustment = "before adjustment"),
                         pivot_longer(as_tibble(ppMat.adj, rownames = "id"), -id, names_to = "sample", values_to = "count") %>% mutate(adjustment = "after adjustment")) %>%
    mutate(adjustment = factor(adjustment, levels = c("before adjustment","after adjustment")),
           count = log2(count)) %>%
    filter(!is.na(count))
  
  # For phosphorylation intensity box plots
  medTab <- group_by(ppPlotTab, adjustment) %>%
    summarise(medVal = median(count, na.rm=TRUE))
  ppBoxplot <- ggplot(ppPlotTab, aes(x=sample, y=count, fill = adjustment)) +
    geom_boxplot() +
    geom_hline(data = medTab, aes(yintercept = medVal), linetype = "dashed", color = "red") +
    facet_wrap(~adjustment, ncol=1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "none") +
    xlab("") + ylab("log2(intensity)") +
    ggtitle("Box plot of phosphorylation in PP samples")
  
  return(list(ratioTrendPlot = ratioTrendPlot,
              ratioBoxplot = ratioBoxplot,
              ppBoxplot = ppBoxplot))
}