#' @name performDifferentialExp
#'
#' @title Perform Differential Expression Analysis
#'
#' @description
#' \code{performDifferentialExp} performs differential expression analysis on a given \code{SummarizedExperiment} object using either the \code{limma} or \code{ProDA} method.
#'
#' @param se A \code{SummarizedExperiment} object containing the data.
#' @param assay A \code{character} string specifying the assay to use for the analysis.
#' @param method A \code{character} string specifying the method to use for differential expression analysis ("limma" or "ProDA"). Default is "limma".
#' @param condition A \code{character} string specifying the condition column in colData(se). Default is \code{NULL}.
#' @param reference A \code{character} string or vector specifying the reference group.
#' @param target A \code{character} string or vector specifying the target group.
#' @param refTime A \code{character} string or vector specifying the reference time points. Default is \code{NULL}.
#' @param targetTime A \code{character} string or vector specifying the target time points. Default is \code{NULL}.
#' @param pairedTtest A \code{logical} value specifying to perform paired t-test or not. Default is \code{FALSE}.
#'
#' @return A list containing:
#' \item{resDE}{A \code{tibble} with the differential expression results.}
#' \item{seSub}{A \code{SummarizedExperiment} object subset to the samples used in the analysis.}
#'
#' @details
#' This function is designed to facilitate differential expression analysis on a \code{SummarizedExperiment} (SE) object. The function allows users to specify various parameters to tailor the analysis to their specific experimental setup.
#'
#' The main steps of the function are as follows:
#'
#' 1. Sample Selection: Based on the provided \code{condition}, \code{reference}, and \code{target} arguments, the function identifies the relevant samples for the analysis. If time points (\code{refTime} and \code{targetTime}) are provided, it further refines the sample selection.
#'
#' 2. Subsetting the SE Object: The SE object is subsetted to include only the selected samples. A new column \code{comparison} is added to the \code{colData}, indicating whether each sample belongs to the reference or target group.
#'
#' 3. Design Matrix Construction: The function constructs a design matrix for the differential expression analysis. If the SE object contains a \code{subjectID} column, this is included in the design to account for repeated measures or paired samples.
#'
#' 4. Differential Expression Analysis: Depending on the specified \code{method}, the function performs the differential expression analysis using either the \code{limma} or \code{ProDA} package:
#'     - \code{Limma}: The function fits a linear model to the expression data and applies empirical Bayes moderation to the standard errors. The results are then extracted and formatted.
#'     - \code{ProDA}: The function fits a probabilistic dropout model to the expression data and tests for differential expression. The results are then extracted and formatted.
#'
#' 5. Result Formatting: The differential expression results are merged with the metadata from the SE object, and the resulting table is formatted into a tibble. The table includes columns for log2 fold change (log2FC), test statistic (stat), p-value (pvalue), adjusted p-value (padj), and gene/feature ID (ID).
#'
#' The function returns a \code{list} containing the formatted differential expression results and the subsetted SE object. This allows users to further explore or visualize the results as needed.
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom proDA proDA test_diff
#' @importFrom dplyr rename select filter arrange as_tibble
#' @importFrom tibble as_tibble
#' @importFrom SummarizedExperiment assays colData
#' @importFrom stats model.matrix
#'
#' @examples
#' library(SummarizedExperiment)
#' # Load multiAssayExperiment object
#' data("dda_example")
#' # Get SummarizedExperiment object
#' se <- dda_example[["Proteome"]]
#' colData(se) <- colData(dda_example)
#' # Preprocess the proteome assay
#' result <- preprocessProteome(se, normalize = TRUE)
#' # Call the function to perform differential expression analyis
#' performDifferentialExp(se = result, assay = "Intensity", method = "limma", reference = "1stCrtl", target = "EGF", condition = "treatment")
#'
#' @export
performDifferentialExp <- function(se, assay, method = "limma", condition = NULL, reference, target, refTime = NULL, targetTime = NULL, pairedTtest = FALSE) {
  # Stop if method is other than limma or proDA
  if (!(method %in% c("limma", "ProDA"))) stop("Invalid method!! Provide either limma or ProDA")

  # Stop if method is other than limma or proDA
  if (!is.null(condition)) {
    if (is.null(se[[condition]])) stop("Invalid condtion! Provide condition which is part of the given Summarized Experiment object")
  }

  # Identify samples based on condition and time points if provided
  if (!is.null(condition)) {
    if (!is.null(refTime)) {
      referenceID <- se[, se[[condition]] %in% reference & se$timepoint %in% refTime]$sample
      targetID <- se[, se[[condition]] %in% target & se$timepoint %in% targetTime]$sample
    } else {
      referenceID <- se[, se[[condition]] %in% reference]$sample
      targetID <- se[, se[[condition]] %in% target]$sample
    }
  } else {
    referenceID <- reference
    targetID <- target
  }

  # Subset the SummarizedExperiment object to the selected samples
  seSub <- se[, se$sample %in% c(referenceID, targetID)]
  # Create a comparison column indicating reference and target samples
  seSub$comparison <- ifelse(seSub$sample %in% referenceID, "reference", "target")
  seSub$comparison <- factor(seSub$comparison, levels = c("reference", "target"))

  # Extract the expression matrix and colData
  seqMat <- seSub
  exprMat <- assays(seqMat)[[assay]]
  colData <- data.frame(colData(seqMat))

  # Create the design matrix for the differential expression analysis
  if (!is.null(seSub$subjectID) && pairedTtest) {
    design <- model.matrix(~ subjectID + comparison, data = colData)
  } else {
    design <- model.matrix(~comparison, data = colData)
  }
  # Get the names of the design matrix columns
  resNames <- colnames(design)
  # Extract the metadata from the SummarizedExperiment object
  meta <- as.data.frame(elementMetadata(seSub))

  # Perform differential expression analysis using the specified method
  if (method == "limma") {
    fit <- limma::lmFit(exprMat, design = design)
    fit2 <- eBayes(fit)
    resDE <- topTable(fit2, number = Inf, coef = resNames[length(resNames)])
    # Merge the results with metadata and format the results
    resDE <- merge(resDE, meta, by = 0, all = TRUE)
    resDE <- as_tibble(resDE) %>%
      dplyr::rename(
        log2FC = logFC, stat = t,
        pvalue = P.Value, padj = adj.P.Val, ID = Row.names
      ) %>%
      select(-c(B, AveExpr)) %>%
      filter(!is.na(padj)) %>%
      arrange(pvalue)
  } else if (method == "ProDA") {
    fit <- proDA::proDA(exprMat, design = design)
    resDE <- proDA::test_diff(fit, contrast = resNames[length(resNames)])
    # Merge the results with metadata and format the results
    rownames(resDE) <- resDE[, 1]
    resDE[, 1] <- NULL
    resDE <- merge(resDE, meta, by = 0, all = TRUE)
    resDE <- as_tibble(resDE) %>%
      dplyr::rename(
        log2FC = diff, stat = t_statistic,
        pvalue = pval, padj = adj_pval, ID = Row.names
      ) %>%
      select(-c(se, df, avg_abundance, n_approx, n_obs)) %>%
      filter(!is.na(padj)) %>%
      arrange(pvalue)
  }

  # Return the results and the subsetted SummarizedExperiment object
  return(list(resDE = resDE, seSub = seSub))
}


#' @name plotVolcano
#'
#' @title Plot Volcano Plot for Differential Expression Analysis
#'
#' @description
#' \code{plotVolcano} generates a volcano plot to visualize differential expression results.
#'
#' @param tableDE A \code{data frame} containing differential expression results with columns 'ID', 'log2FC', 'pvalue', and 'Gene'.
#' @param pFilter A \code{numeric} value specifying the p-value threshold for significance. Default is 0.05.
#' @param fcFilter A \code{numeric} value specifying the log2 fold-change threshold for significance. Default is 0.5.
#'
#' @return A \code{ggplot2} object representing the volcano plot.
#'
#' @details
#' This function creates a volcano plot where differentially expressed genes are categorized as 'Up', 'Down', or 'Not Sig' based on the provided p-value and log2 fold-change thresholds. Points on the plot are color-coded to indicate their expression status.
#'
#' @importFrom dplyr mutate case_when
#' @import ggplot2
#'
#' @examples
#' library(SummarizedExperiment)
#' # Load multiAssayExperiment object
#' data("dda_example")
#' # Get SummarizedExperiment object
#' se <- dda_example[["Proteome"]]
#' colData(se) <- colData(dda_example)
#' # Preprocess the proteome assay
#' result <- preprocessProteome(se, normalize = TRUE)
#' # Call the function to perform differential expression analyis
#' de <- performDifferentialExp(se = result, assay = "Intensity", method = "limma", reference = "1stCrtl", target = "EGF", condition = "treatment")
#' # Plot the volcano plot from the result
#' plotVolcano(de$resDE)
#'
#' @export
plotVolcano <- function(tableDE, pFilter = 0.05, fcFilter = 0.5) {
  if (!("log2FC" %in% colnames(tableDE))) stop("column 'log2FC' not found")
  if (!("pvalue" %in% colnames(tableDE))) stop("column 'pvalue' not found")
  if (!("Gene" %in% colnames(tableDE))) stop("column 'Gene' not found")
  if (!("ID" %in% colnames(tableDE))) stop("column 'ID' not found")

  # Convert the input table to a data frame and ensure the 'ID' column is of type character
  dataVolcano <- data.frame(tableDE)
  dataVolcano$ID <- as.character(dataVolcano$ID)
  # Categorize each gene based on the provided p-value and log2 fold-change thresholds
  dataVolcano <- mutate(dataVolcano, expression = case_when(
    dataVolcano$log2FC >= as.numeric(fcFilter) & dataVolcano$pvalue <= as.numeric(pFilter) ~ "Up",
    dataVolcano$log2FC <= -as.numeric(fcFilter) & dataVolcano$pvalue <= as.numeric(pFilter) ~ "Down",
    dataVolcano$pvalue > as.numeric(pFilter) | (dataVolcano$log2FC < as.numeric(fcFilter) & dataVolcano$log2FC > -as.numeric(fcFilter)) ~ "Not Sig"
  ))

  # Create the volcano plot
  v <- ggplot(dataVolcano, aes(x = log2FC, y = -log10(pvalue), label = Gene, customdata = ID)) +
    # Add vertical lines for fold-change thresholds
    geom_vline(xintercept = 0, color = "black", linetype = "solid", linewidth = 0.25) +
    geom_vline(xintercept = as.numeric(fcFilter), color = "darkgrey", linetype = "dashed") +
    geom_vline(xintercept = -as.numeric(fcFilter), color = "darkgrey", linetype = "dashed") +
    # Add horizontal lines for p-value thresholds
    geom_hline(yintercept = -log10(as.numeric(pFilter)), color = "darkgrey", linetype = "dashed") +
    annotate(
      x = 5.0, y = -log10(as.numeric(pFilter)) - 0.1, label = paste("P-value = ", as.numeric(pFilter)),
      geom = "text", size = 3, color = "darkgrey"
    ) +
    geom_hline(yintercept = -log10(0.25), color = "darkgrey", linetype = "dashed") +
    annotate(
      x = 5.0, y = 0.5, label = paste("P-value = ", 0.25),
      geom = "text", size = 3, color = "darkgrey"
    ) +
    # Plot the points and color them based on their expression status
    geom_point(aes(color = expression), size = 0.9) +
    scale_color_manual(values = c("Up" = "firebrick3", "Down" = "navy", "Not Sig" = "darkgrey")) +
    xlab("absolute log2(Quantity) difference") +
    ggtitle("Volcano plot") +
    theme(plot.title = element_text(hjust = 0.5))
  return(v)
}

#' @name intensityBoxPlot
#'
#' @title Plot Boxplot of Intensity Data
#'
#' @description
#' \code{intensityBoxPlot} creates a boxplot for the Intensity data of a given gene or feature, with optional subject-specific lines.
#'
#' @param se A \code{SummarizedExperiment} object containing the data.
#' @param id \code{Character}. The identifier of the gene or feature to plot.
#' @param symbol \code{Character}. The symbol or name of the gene or feature to use as the plot title.
#'
#' @return A \code{ggplot2} object representing the boxplot of the intensity data.
#'
#' @details
#' This function generates a boxplot for the intensity data of a specified gene or feature from a \code{SummarizedExperiment} (SE) object. The plot shows the distribution of normalized intensities across different groups specified in the \code{comparison} column of the SE object.
#'
#' The function can handle both grouped data and repeated measures:
#' - If the SE object does not contain a \code{subjectID} column, the function plots a standard boxplot grouped by the \code{comparison} column.
#' - If the SE object contains a \code{subjectID} column, the function adds lines connecting the points for each subject across the groups, providing a visual indication of subject-specific changes.
#'
#' The \code{boxplot} is customized with various aesthetic elements, such as box width, transparency, point size, axis labels, and title formatting.
#'
#' @importFrom SummarizedExperiment assays
#' @import ggplot2
#'
#' @examples
#' library(SummarizedExperiment)
#' # Load multiAssayExperiment object
#' data("dda_example")
#' # Get SummarizedExperiment object
#' se <- dda_example[["Proteome"]]
#' colData(se) <- colData(dda_example)
#' # Preprocess the proteome assay
#' result <- preprocessProteome(se, normalize = TRUE)
#' # Call the function to perform differential expression analyis
#' de <- performDifferentialExp(se = result, assay = "Intensity", method = "limma", reference = "1stCrtl", target = "EGF", condition = "treatment")
#' # Plot the box plot for the given id and symbol
#' intensityBoxPlot(de$seSub, "p99", "PPP6C")
#'
#' @export
intensityBoxPlot <- function(se, id, symbol) {
  exprMat <- assays(se)[["Intensity"]]

  # Check if the SE object contains subject-specific data
  if (is.null(se$subjectID)) {
    # Prepare data frame for plotting without subject-specific information
    plotTab <- data.frame(
      group = se$comparison,
      value = exprMat[id, ]
    )
    p <- ggplot(plotTab, aes(x = group, y = value))
  } else {
    # Prepare data frame for plotting with subject-specific information
    plotTab <- data.frame(
      group = se$comparison,
      value = exprMat[id, ],
      subjectID = se$subjectID
    )
    p <- ggplot(plotTab, aes(x = group, y = value)) +
      geom_line(aes(group = subjectID), linetype = "dotted", color = "grey50")
  }

  # Create the boxplot with additional formatting
  p <- p + geom_boxplot(aes(fill = group),
    width = 0.5, alpha = 0.5,
    outlier.shape = NA
  ) +
    geom_point() +
    ylab("Normalized Intensities") + xlab("") +
    ggtitle(symbol) + theme_bw() +
    theme(
      text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 15)
    )

  return(p)
}
