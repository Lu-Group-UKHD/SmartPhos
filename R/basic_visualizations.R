#' @name plotMissing
#'
#' @title Plot Missing Data Completeness
#'
#' @description
#' \code{plotMissing} generates a bar plot showing the completeness (percentage
#' of non-missing values) for each sample in a \code{SummarizedExperiment}
#' object.
#'
#' @param se A \code{SummarizedExperiment} object containing the assay data.
#'
#' @return A \code{ggplot2} object showing the percentage of completeness for
#' each sample.
#'
#' @details
#' This function calculates the percentage of non-missing values for each sample
#' in the provided \code{SummarizedExperiment} object. It then generates a bar
#' plot where each bar represents a sample, and the height of the bar
#' corresponds to the completeness (percentage of non-missing values) of that
#' sample.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble tibble
#' @import ggplot2
#'
#' @examples
#' library(SummarizedExperiment)
#' # Load multiAssayExperiment object
#' data("dda_example")
#' # Get SummarizedExperiment object
#' se <- dda_example[["Phosphoproteome"]]
#' colData(se) <- colData(dda_example)
#' # Call the function
#' plotMissing(se)
#'
#' @export
plotMissing <- function(se) {

  perNA <- NULL
  # Extract the assay data from the SummarizedExperiment object
  countMat <- assay(se)

  # Create a table with sample names and their corresponding percentage of
  # non-missing values
  plotTab <- tibble(
    sample = se$sample,
    perNA = colSums(is.na(countMat)) / nrow(countMat)
  )

  # Generate the bar plot using ggplot2
  missPlot <- ggplot(plotTab, aes(x = sample, y = 1 - perNA)) +
    geom_bar(stat = "identity") +
    ggtitle("Percentage of sample completeness") +
    ylab("completeness") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  return(missPlot)
}


#' @name plotIntensity
#'
#' @title Plot Intensity Boxplots
#'
#' @description
#' \code{plotIntensity} generates boxplots of assay intensities for each sample
#' in a \code{SummarizedExperiment} object. Optionally, the boxplots can be
#' colored based on a specified metadata column. The function handles missing
#' values by filtering them out before plotting.
#'
#' @param se A \code{SummarizedExperiment} object containing the assay data and
#' metadata.
#' @param colorByCol A \code{character} string specifying the metadata column
#' to use for coloring the boxplots. Default is "none".
#'
#' @return A \code{ggplot2} object showing boxplots of intensities for each
#' sample.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr filter left_join
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @import ggplot2
#'
#' @examples
#' library(SummarizedExperiment)
#' # Load multiAssayExperiment object
#' data("dia_example")
#' # Get SummarizedExperiment object
#' se <- dia_example[["Phosphoproteome"]]
#' colData(se) <- colData(dia_example)
#' # Preprocess the phosphoproteome assay
#' result <- preprocessPhos(seData = se, normalize = TRUE, impute = "QRILC")
#' # Call the plotting function
#' plotIntensity(result, colorByCol = "replicate")
#'
#' @export
plotIntensity <- function(se, colorByCol = "none") {

  id <- value <- name <- NULL

  # Extract the assay data from the SummarizedExperiment object
  countMat <- assay(se)

  # Convert assay data to a tibble, pivot to long format, and
  # filter out missing values
  countTab <- countMat %>%
    as_tibble(rownames = "id") %>%
    pivot_longer(-id) %>%
    filter(!is.na(value))

  # Extract metadata from the SummarizedExperiment object
  meta <- as.data.frame(colData(se))
  # Join the count data with metadata
  countTabmeta <- left_join(countTab, meta, by = c("name" = "sample"))

  # Create the ggplot object with boxplots of intensities
  g <- ggplot(countTabmeta, aes(x = name, y = value)) +
    ggtitle("Boxplot of intensities") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )

  # Add color to the boxplots if a valid metadata column is specified
  if (colorByCol == "none") {
    g <- g + geom_boxplot()
  } else {
    g <- g + geom_boxplot(aes(color = !!sym(colorByCol)))
  }

  return(g)
}

#' @name plotPCA
#'
#' @title Plot PCA
#'
#' @description
#' \code{plotPCA} generates a PCA plot using the results from a PCA analysis
#' and a \code{SummarizedExperiment} object. The points on the plot can be
#' colored and shaped based on metadata.
#'
#' @param pca A PCA result object, typically obtained from \code{prcomp}.
#' @param se A \code{SummarizedExperiment} object containing the metadata.
#' @param xaxis A \code{character} string specifying which principal component
#' to use for the x-axis. Default is "PC1".
#' @param yaxis A \code{character} string specifying which principal component
#' to use for the y-axis. Default is "PC2".
#' @param color A \code{character} string specifying the metadata column to use
#' for coloring the points. Default is "none".
#' @param shape A \code{character} string specifying the metadata column to use
#' for shaping the points. Default is "none".
#'
#' @return A \code{ggplot2} object showing the PCA plot.
#'
#' @details
#' This function creates a PCA plot using the scores from a PCA result object
#' and metadata from a \code{SummarizedExperiment} object. The x-axis and y-axis
#' can be customized to display different principal components, and the points
#' can be optionally colored and shaped based on specified metadata columns.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column
#' @importFrom rlang sym
#' @import ggplot2
#'
#' @examples
#' # Load multiAssayExperiment object
#' data("dia_example")
#' # Get SummarizedExperiment object
#' se <- dia_example[["Phosphoproteome"]]
#' SummarizedExperiment::colData(se) <- SummarizedExperiment::colData(
#' dia_example)
#' # Generate the imputed assay
#' result <- preprocessPhos(seData = se, normalize = TRUE, impute = "QRILC")
#' # Perform PCA
#' pcaResult <- stats::prcomp(t(
#' SummarizedExperiment::assays(result)[["imputed"]]),
#' center = TRUE, scale. = TRUE)
#' # Plot PCA results
#' plotPCA(pca = pcaResult, se = result, color = "treatment")
#'
#' @export
plotPCA <- function(pca, se, xaxis = "PC1", yaxis = "PC2", color = "none",
                    shape = "none") {
  # Calculate the proportion of variance explained by each principal component
  varExplained <- pca$sdev^2 / sum(pca$sdev^2)
  # Convert the PCA result to a data frame
  pcaDf <- as.data.frame(pca[["x"]])
  # Convert the metadata to a data frame
  meta <- as.data.frame(colData(se))
  # Join the PCA scores with the metadata
  pcaMeta <- left_join(rownames_to_column(pcaDf),
    meta,
    by = c("rowname" = "sample")
  )

  # Create the initial ggplot object with labels for variance explained
  g <- ggplot(pcaMeta, aes(
    x = !!sym(xaxis), y = !!sym(yaxis),
    text = paste("sample:", meta$sample)
  )) +
    theme_bw() +
    theme(legend.position = "top") +
    labs(
      x = paste0(
        xaxis, ": ",
        round(varExplained[as.numeric(strsplit(xaxis, "PC")[[1]][2])] * 100,
              1), "%"),
      y = paste0(
        yaxis, ": ",
        round(varExplained[as.numeric(strsplit(yaxis, "PC")[[1]][2])] * 100,
              1), "%")
    ) +
    scale_shape(solid = FALSE)

  # Add points to the plot with optional color and shape aesthetics
  if (color == "none" & shape == "none") {
    g <- g + geom_point(size = 2)
  } else if (color == "none") {
    g <- g + geom_point(aes(shape = !!sym(shape)), size = 2)
  } else if (shape == "none") {
    g <- g + geom_point(aes(color = !!sym(color)), size = 2)
  } else {
    g <- g + geom_point(
      aes(
        color = !!sym(color),
        shape = !!sym(shape)
      ),
      size = 2
    )
  }
}

#' @name plotHeatmap
#'
#' @title Plot Heatmap of Intensity assay
#'
#' @description
#' \code{plotHeatmap} generates a heatmap for intensity assay for different
#' conditions, including top variants, differentially expressed genes, and
#' selected time series clusters.
#'
#' @param type A \code{character} string specifying the type of heatmap to plot.
#' Options are "Top variant", "Differentially expressed", and "Selected time
#' series cluster".
#' @param se A \code{SummarizedExperiment} object containing the imputed
#' intensity assay.
#' @param data An optional \code{data frame} containing additional data for
#' "Differentially expressed" and "Selected time series cluster" types. Default
#' is \code{NULL}.
#' @param top A \code{numeric} value specifying the number of top variants to
#' plot. Default is 100.
#' @param cutCol A \code{numeric} value specifying the number of clusters for
#' columns. Default is 1.
#' @param cutRow A \code{numeric} value specifying the number of clusters for
#' rows. Default is 1.
#' @param clustCol A \code{logical} value indicating whether to cluster columns.
#' Default is \code{TRUE}.
#' @param clustRow A \code{logical} value indicating whether to cluster rows.
#' Default is \code{TRUE}.
#' @param annotationCol A \code{character} vector specifying the columns in the
#' metadata to use for annotation. Default is \code{NULL}.
#' @param title A \code{character} string specifying the title of the heatmap.
#' Default is \code{NULL}.
#'
#' @return A \code{pheatmap} object showing the heatmap of Intensity data.
#'
#' @details
#' This function creates a heatmap using the Intensity assay from a
#' \code{SummarizedExperiment} object. The heatmap can show the top variants
#' based on standard deviation, differentially expressed genes, or selected time
#' series clusters. Row normalization is performed, and the heatmap can include
#' annotations based on specified metadata columns.
#'
#' @importFrom SummarizedExperiment assays colData rowData
#' @importFrom dplyr arrange left_join
#' @importFrom pheatmap pheatmap
#' @importFrom tibble rownames_to_column
#' @importFrom grDevices colorRampPalette
#' @importFrom rlang sym
#' @import ggplot2
#'
#' @examples
#' library(SummarizedExperiment)
#' # Load multiAssayExperiment object
#' data("dia_example")
#' # Get SummarizedExperiment object
#' se <- dia_example[["Phosphoproteome"]]
#' colData(se) <- colData(dia_example)
#' # Generate the imputed assay
#' result <- preprocessPhos(seData = se, normalize = TRUE, impute = "QRILC")
#' # Plot heatmap for top variant
#' plotHeatmap(type = "Top variant", top = 10, se = result, cutCol = 2)
#'
#' @export
plotHeatmap <- function(type = c("Top variant", "Differentially expressed",
                                 "Selected time series cluster"),
                        se, data = NULL, top = 100, cutCol = 1,
                        cutRow = 1, clustCol = TRUE, clustRow = TRUE,
                        annotationCol = NULL, title = NULL) {

    type <- match.arg(type)
    # Select the appropriate intensity assay and gene IDs based on the type of
    # heatmap
    if (type == "Top variant") {
        exprMat <- assays(se)[["imputed"]]
        sds <- apply(exprMat, 1, sd)
        orderID <- names(sort(sds, decreasing = TRUE))
        geneIDs <- orderID[seq(1, as.integer(top))]
        exprMat <- exprMat[geneIDs, ]
        geneSymbol <- rowData(se[match(geneIDs, rownames(se)), ])$Gene
    } else if (type == "Differentially expressed") {
        if (!is.null(data)) {
            geneIDs <- arrange(data, stat)$ID
            exprMat <- assays(se)[["imputed"]][geneIDs, ]
            geneSymbol <- data[match(geneIDs, data$ID), ]$Gene
        } else {
            message("Please give data argument")
        }
    } else if (type == "Selected time series cluster") {
        if (!is.null(data)) {
            geneIDs <- unique(data$ID)
            exprMat <- assays(se)[["imputed"]][geneIDs, ]
            geneSymbol <- data[match(geneIDs, data$ID), ]$Gene
        } else {
            message("Please give data argument")
        }
    }

    # Prepare column annotations from the metadata
    if (!is.null(annotationCol)) {
        cd <- as.data.frame(colData(se))
        annCol <- cd[row.names(cd) %in% colnames(exprMat), ][c(annotationCol)]
        row.names(annCol) <- colnames(exprMat)
    } else {
        annCol <- NA
    }


    # Prepare the title of heatmap if Null
    if (is.null(title))
        title <- type

    # Prepare color scale for the heatmap
    color <- colorRampPalette(c("navy", "white", "firebrick"))(100)

    # Perform row normalization and clip extreme values
    exprMat <- t(scale(t(exprMat)))
    exprMat[exprMat > 4] <- 4
    exprMat[exprMat < -4] <- -4

    # Plot the heatmap based on the type and whether annotations are provided
    if (type == "Top variant") {
        p <- pheatmap(exprMat, color = color, labels_row = geneSymbol,
                      treeheight_row = 0, treeheight_col = 0, main = title,
                      cutree_cols = cutCol, cutree_rows = cutRow,
                      annotation_col = annCol)
    } else {
        # Sort the columns by their names before plotting
        p <- pheatmap(exprMat, color = color, labels_row = geneSymbol,
                      treeheight_row = 0, treeheight_col = 0, main = title,
                      cluster_rows = clustRow, cluster_cols = clustCol,
                      cutree_cols = cutCol, cutree_rows = cutRow,
                      annotation_col = annCol)
    }
    return(p)
}
