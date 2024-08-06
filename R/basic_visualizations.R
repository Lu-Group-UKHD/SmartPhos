#' @name plotMissing 
#' 
#' @title Plot Missing Data Completeness
#'
#' @description 
#' `plotMissing` generates a bar plot showing the completeness (percentage of non-missing values) for each sample in a SummarizedExperiment object.
#'
#' @param se A SummarizedExperiment object containing the assay data.
#'
#' @return A ggplot object showing the percentage of completeness for each sample.
#'
#' @details
#' This function calculates the percentage of non-missing values for each sample in the provided SummarizedExperiment object. It then generates a bar plot where each bar represents a sample, and the height of the bar corresponds to the completeness (percentage of non-missing values) of that sample.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom ggplot2 ggplot aes geom_bar ggtitle ylab theme element_text
#' @importFrom tibble tibble
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object:
#' plot <- plotMissing(se)
#' print(plot)
#'
#' @export
plotMissing <- function(se) {
  
  # Extract the assay data from the SummarizedExperiment object
  countMat <- assay(se)
  
  # Create a table with sample names and their corresponding percentage of non-missing values
  plotTab <- tibble(sample = se$sample, 
                    perNA = colSums(is.na(countMat))/nrow(countMat))
  
  # Generate the bar plot using ggplot2
  missPlot <- ggplot(plotTab, aes(x = sample, y = 1-perNA)) +
    geom_bar(stat = "identity") +
    ggtitle("Percentage of sample completeness") +
    ylab("completeness") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  return(missPlot)
}


#' @name plotIntensity
#' 
#' @title Plot Intensity Boxplots
#'
#' @description
#' `plotIntensity` generates boxplots of assay intensities for each sample in a SummarizedExperiment object. Optionally, the boxplots can be colored based on a specified metadata column. The function handles missing values by filtering them out before plotting.
#'
#' @param se A SummarizedExperiment object containing the assay data and metadata.
#' @param color A character string specifying the metadata column to use for coloring the boxplots. Default is "none".
#'
#' @return A ggplot object showing boxplots of intensities for each sample.
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom ggplot2 ggplot aes geom_boxplot ggtitle theme element_text
#' @importFrom dplyr filter left_join
#' @importFrom tidyr pivot_longer
#' @importFrom tibble as_tibble
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object:
#' plot <- plotIntensity(se, color = "group")
#' print(plot)
#'
#' @export
plotIntensity <- function(se, color = "none") {
  
  # Extract the assay data from the SummarizedExperiment object
  countMat <- assay(se)
  
  # Convert the assay data to a tibble, pivot to long format, and filter out missing values
  countTab <- countMat %>% as_tibble(rownames = "id") %>% 
    pivot_longer(-id) %>%
    filter(!is.na(value))
  
  # Extract metadata from the SummarizedExperiment object
  meta <- as.data.frame(colData(se))
  # Join the count data with metadata
  countTabmeta <- left_join(countTab, meta, by = c('name' = 'sample'))
  
  # Create the ggplot object with boxplots of intensities
  g <- ggplot(countTabmeta, aes(x = name, y = value)) +
    ggtitle("Boxplot of intensities") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0),
          plot.title = element_text(hjust = 0.5, face = "bold"))
  
  # Add color to the boxplots if a valid metadata column is specified
  if (color == "none"){
    g <- g + geom_boxplot()
  }
  else {
    g <- g + geom_boxplot(aes_string(fill = color))
  }
  
  return(g)
}

#' @name plotPCA
#' 
#' @title Plot PCA
#'
#' @description
#' `plotPCA` generates a PCA plot using the results from a PCA analysis and a SummarizedExperiment object. The points on the plot can be colored and shaped based on metadata.
#'
#' @param pca A PCA result object, typically obtained from \code{prcomp}.
#' @param se A SummarizedExperiment object containing the metadata.
#' @param xaxis A character string specifying which principal component to use for the x-axis. Default is `"PC1"`.
#' @param yaxis A character string specifying which principal component to use for the y-axis. Default is `"PC2"`.
#' @param color A character string specifying the metadata column to use for coloring the points. Default is `"none"`.
#' @param shape A character string specifying the metadata column to use for shaping the points. Default is `"none"`.
#'
#' @return A ggplot object showing the PCA plot.
#'
#' @details
#' This function creates a PCA plot using the scores from a PCA result object and metadata from a SummarizedExperiment object. The x-axis and y-axis can be customized to display different principal components, and the points can be optionally colored and shaped based on specified metadata columns.
#'
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes geom_point theme_bw theme labs scale_shape
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column
#' @importFrom rlang sym
#' @examples
#' # Assuming 'pca' is a PCA result object and 'se' is a SummarizedExperiment object:
#' plot <- plotPCA(pca, se, xaxis = "PC1", yaxis = "PC2", color = "group", shape = "type")
#' print(plot)
#'
#' @export
plotPCA <- function(pca, se, xaxis = "PC1", yaxis = "PC2", color = "none", shape = "none") {
  
  # Calculate the proportion of variance explained by each principal component
  varExplained <- pca$sdev^2/sum(pca$sdev^2)
  # Convert the PCA result to a data frame
  pcaDf <- as.data.frame(pca[["x"]])
  # Convert the metadata to a data frame
  meta <- as.data.frame(colData(se))
  # Join the PCA scores with the metadata 
  pcaMeta <- left_join(rownames_to_column(pcaDf),
                       meta, by = c("rowname" = "sample"))
  
  # Create the initial ggplot object with labels for variance explained
  g <- ggplot(pcaMeta, aes(x = !!sym(xaxis), y = !!sym(yaxis),
                           text = paste("sample:", meta$sample))) +
    theme_bw() +
    theme(legend.position="top") +
    labs(x=paste0(xaxis,": ",
                  round(varExplained[as.numeric(strsplit(xaxis, "PC")[[1]][2])]*100, 1), "%"),
         y=paste0(yaxis,": ",
                  round(varExplained[as.numeric(strsplit(yaxis, "PC")[[1]][2])]*100, 1), "%")) +
    scale_shape(solid = FALSE)
  
  # Add points to the plot with optional color and shape aesthetics
  if (color == "none" & shape == "none") {
    g <- g + geom_point(size = 2)
  }
  else if (color == "none") {
    g <- g + geom_point(aes(shape = !!sym(shape)), size = 2)
  }
  else if (shape == "none") {
    g <- g + geom_point(aes(color = !!sym(color)), size = 2)
  }
  else {
    g <- g + geom_point(aes(color = !!sym(color),
                                   shape = !!sym(shape)),
                        size = 2)
  }
}

#' @name plotHeatmap
#' 
#' @title Plot Heatmap of Intensity assay
#'
#' @description
#' `plotHeatmap` generates a heatmap for intensity assay for different conditions, including top variants, differentially expressed genes, and selected time series clusters.
#'
#' @param type A character string specifying the type of heatmap to plot. Options are `"Top variant"`, `"Differentially expressed"`, and `"Selected time series cluster"`.
#' @param se A SummarizedExperiment object containing the imputed intensity assay.
#' @param data An optional data frame containing additional data for `"Differentially expressed"` and `"Selected time series cluster"` types. Default is `NULL`.
#' @param top An integer specifying the number of top variants to plot. Default is `100`.
#' @param cutCol An integer specifying the number of clusters for columns. Default is `1`.
#' @param cutRow An integer specifying the number of clusters for rows. Default is `1`.
#' @param clustCol A logical value indicating whether to cluster columns. Default is `TRUE`.
#' @param clustRow A logical value indicating whether to cluster rows. Default is `TRUE`.
#' @param annotationCol A character vector specifying the columns in the metadata to use for annotation. Default is `NULL`.
#' @param title A character string specifying the title of the heatmap.
#'
#' @return A pheatmap object showing the heatmap of Intensity data.
#'
#' @details
#' This function creates a heatmap using the Intensity assay from a SummarizedExperiment object. The heatmap can show the top variants based on standard deviation, differentially expressed genes, or selected time series clusters. Row normalization is performed, and the heatmap can include annotations based on specified metadata columns.
#'
#' @importFrom SummarizedExperiment assays colData rowData
#' @importFrom ggplot2 ggplot aes geom_point theme_bw theme labs scale_shape
#' @importFrom dplyr arrange left_join
#' @importFrom pheatmap pheatmap
#' @importFrom tibble rownames_to_column
#' @importFrom rlang sym
#' @examples
#' # Assuming 'se' is a SummarizedExperiment object and 'data' is a data frame:
#' heatmap <- plotHeatmap("Top variant", se, top = 50, annotationCol = "group", title = "Top Variants Heatmap")
#' print(heatmap)
#'
#' @export
plotHeatmap <- function(type, se, data = NULL, top = 100, cutCol = 1, cutRow = 1, clustCol = TRUE, clustRow = TRUE, annotationCol, title) {
  
  # Select the appropriate intensity assay and gene IDs based on the type of heatmap
  if (type == "Top variant") {
    exprMat <- assays(se)[["imputed"]]
    sds <- apply(exprMat, 1, sd)
    orderID <- names(sort(sds, decreasing = TRUE))
    geneIDs <- orderID[seq(1, as.integer(top))]
    exprMat <- exprMat[geneIDs,]
    geneSymbol <- rowData(se[match(geneIDs, rownames(se)),])$Gene
  } 
  else if (type == "Differentially expressed") {
    if(!is.null(data)) {
      geneIDs <- arrange(data, stat)$ID
      exprMat <- assays(se)[["imputed"]][geneIDs,] 
      geneSymbol <- data[match(geneIDs, data$ID),]$Gene
    }
    else {
      print("Please give data argument")
    }
  }
  else if (type == "Selected time series cluster") {
    if(!is.null(data)) {
      geneIDs <- unique(data$ID)
      exprMat <- assays(se)[["imputed"]][geneIDs,] 
      geneSymbol <- data[match(geneIDs, data$ID),]$Gene
    }
    else {
      print("Please give data argument")
    }
  }
  
  # Prepare column annotations from the metadata
  cd <- as.data.frame(colData(se))
  annCol <- cd[row.names(cd) %in% colnames(exprMat),][c(annotationCol)]
  row.names(annCol) <- colnames(exprMat)
  
  # Prepare color scale for the heatmap
  color <- colorRampPalette(c("navy", "white", "firebrick"))(100)
  
  # Perform row normalization and clip extreme values
  exprMat <- t(scale(t(exprMat)))
  exprMat[exprMat > 4] <- 4
  exprMat[exprMat < -4] <- -4
  
  # Plot the heatmap based on the type and whether annotations are provided
  if (type == "Top variant") {
    if (is.null(annotationCol)) {
      p <- pheatmap(exprMat, color = color,
                    labels_row = geneSymbol,
                    treeheight_row = 0, treeheight_col = 0,
                    main = title,
                    cutree_cols = cutCol,
                    cutree_rows = cutRow)
    }
    else {
      p <- pheatmap(exprMat, color = color,
                    labels_row = geneSymbol,
                    treeheight_row = 0, treeheight_col = 0,
                    main = title,
                    cutree_cols = cutCol,
                    cutree_rows = cutRow,
                    annotation_col = annCol)
    }
  }
  else {
    # Sort the columns by their names before plotting
    exprMat <- exprMat[, sort(colnames(exprMat))]
    if (is.null(annotationCol)) {
      p <- pheatmap(exprMat, color = color,
                    labels_row = geneSymbol,
                    treeheight_row = 0, treeheight_col = 0,
                    main = title,
                    cluster_rows = clustRow, cluster_cols = clustCol,
                    cutree_cols = cutCol,
                    cutree_rows = cutRow)
    }
    else {
      p <- pheatmap(exprMat, color = color,
                    labels_row = geneSymbol,
                    treeheight_row = 0, treeheight_col = 0,
                    main = title,
                    cluster_rows = clustRow, cluster_cols = clustCol,
                    cutree_cols = cutCol,
                    cutree_rows = cutRow,
                    annotation_col = annCol)
    }
  }
  return(p)
}