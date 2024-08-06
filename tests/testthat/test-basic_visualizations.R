# Helper function to create a mock SummarizedExperiment
create_mock_se <- function() {
  
  # Create assay data
  assay_data <- matrix(rnorm(1000), nrow = 100, ncol = 10)
  colnames(assay_data) <- paste0("Sample", 1:10)
  rownames(assay_data) <- paste0("Gene", 1:100)
  
  # Create sample annotations
  sample_data <- data.frame(sample = colnames(assay_data), group = rep(c("A", "B"), each = 5), type = rep(c("X", "Y"), 5))
  rownames(sample_data) <- colnames(assay_data)
  
  # Create gene metadata
  gene_data <- data.frame(Gene = rownames(assay_data), stat = runif(100), ID = paste0("Gene", 1:100))
  rownames(gene_data) <- rownames(assay_data)
  
  # Create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(imputed = assay_data), colData = sample_data, rowData = gene_data)
  
  return(se)
}


######################## Tests for plotMissing() ###############################

test_that("plotMissing calculates completeness correctly", {
  se <- create_mock_se()

  # Extract the assay data
  assay_data <- assay(se)

  # Calculate expected completeness
  expected_completeness <- colSums(is.na(assay_data)) / nrow(assay_data)

  # Calculate completeness using plotMissing function
  plot <- plotMissing(se)
  plot_data <- plot$data

  expect_equal(plot_data$sample, colnames(assay_data))
  expect_equal(plot_data$perNA, expected_completeness)
})

test_that("plotMissing generates the correct ggplot object", {
  se <- create_mock_se()

  plot <- plotMissing(se)

  # Check if the returned object is a ggplot object
  expect_true(is.ggplot(plot))

})

test_that("plotMissing handles SummarizedExperiment with all missing values", {
  assay_data <- matrix(NA, nrow = 3, ncol = 5)
  colnames(assay_data) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")

  sample_data <- data.frame(sample = colnames(assay_data))
  rownames(sample_data) <- colnames(assay_data)

  se <- SummarizedExperiment(assays = list(counts = assay_data), colData = sample_data)

  # Calculate expected completeness
  expected_completeness <- colSums(is.na(assay_data)) / nrow(assay_data)

  plot <- plotMissing(se)
  plot_data <- plot$data

  # Expect completeness to be 0 for all samples
  expect_equal(plot_data$perNA, expected_completeness)
})


######################## Tests for plotIntensity() #############################

test_that("plotIntensity handles missing values correctly", {
  se <- create_mock_se()

  plot <- plotIntensity(se, color = "none")

  # Extract data used in the plot
  plot_data <- ggplot_build(plot)$data[[1]]

  # Ensure no missing values are present in the plot data
  expect_false(any(is.na(plot_data$y)))
})

test_that("plotIntensity generates the correct ggplot object", {
  se <- create_mock_se()

  plot <- plotIntensity(se, color = "none")

  # Check if the returned object is a ggplot object
  expect_true(is.ggplot(plot))

  # Check plot title
  expect_equal(plot$labels$title, "Boxplot of intensities")

  # Check y-axis label
  expect_equal(plot$labels$y, "value")

  # Check x-axis labels
  expect_equal(plot$labels$x, "name")
})

test_that("plotIntensity applies coloring correctly", {
  se <- create_mock_se()

  plot <- plotIntensity(se, color = "group")

  # Check if the fill aesthetic is set correctly
  expect_true("fill" %in% names(plot$labels))
  expect_equal(plot$labels$fill, "group")
})

test_that("plotIntensity works correctly with color set to 'none'", {
  se <- create_mock_se()

  plot <- plotIntensity(se, color = "none")

  # Check if the fill aesthetic is not set
  expect_false("fill" %in% names(plot$labels))
})

test_that("plotIntensity handles SummarizedExperiment with all non-missing values", {
  assay_data <- matrix(1:15, nrow = 3, ncol = 5)
  colnames(assay_data) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")

  sample_data <- data.frame(sample = colnames(assay_data), group = c("A", "B", "A", "B", "A"))
  rownames(sample_data) <- colnames(assay_data)

  se <- SummarizedExperiment(assays = list(counts = assay_data), colData = sample_data)

  plot <- plotIntensity(se, color = "none")
  plot_data <- ggplot_build(plot)$data[[1]]

  # Expect no missing values in the plot data
  expect_false(any(is.na(plot_data$y)))
})

test_that("plotIntensity handles SummarizedExperiment with all missing values", {
  assay_data <- matrix(NA, nrow = 3, ncol = 5)
  colnames(assay_data) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")

  sample_data <- data.frame(sample = colnames(assay_data), group = c("A", "B", "A", "B", "A"))
  rownames(sample_data) <- colnames(assay_data)

  se <- SummarizedExperiment(assays = list(counts = assay_data), colData = sample_data)

  plot <- plotIntensity(se, color = "none")
  plot_data <- ggplot_build(plot)$data[[1]]

  # Expect no data in the plot since all values are missing
  expect_equal(nrow(plot_data), 0)
})


########################### Tests for plotPCA() ################################

test_that("plotPCA generates a ggplot object", {
  se <- create_mock_se()
  pca <- prcomp(t(assay(se)), scale = TRUE)

  plot <- plotPCA(pca, se)

  # Check if the returned object is a ggplot object
  expect_true(is.ggplot(plot))
})

test_that("plotPCA correctly labels the axes with variance explained", {
  se <- create_mock_se()
  pca <- prcomp(t(assay(se)), scale = TRUE)

  plot <- plotPCA(pca, se, xaxis = "PC1", yaxis = "PC2")

  # Extract axis labels
  x_label <- plot$labels$x
  y_label <- plot$labels$y

  # Calculate expected variance explained
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  expected_x_label <- paste0("PC1: ", round(var_explained[1] * 100, 1), "%")
  expected_y_label <- paste0("PC2: ", round(var_explained[2] * 100, 1), "%")

  # Check if the axis labels match the expected labels
  expect_equal(x_label, expected_x_label)
  expect_equal(y_label, expected_y_label)
})

test_that("plotPCA applies coloring correctly", {
  se <- create_mock_se()
  pca <- prcomp(t(assay(se)), scale = TRUE)

  plot <- plotPCA(pca, se, color = "group")

  # Check if the color aesthetic is set correctly
  expect_true("colour" %in% names(plot$labels))
  expect_equal(plot$labels$colour, "group")
})

test_that("plotPCA applies shaping correctly", {
  se <- create_mock_se()
  pca <- prcomp(t(assay(se)), scale = TRUE)

  plot <- plotPCA(pca, se, shape = "type")

  # Check if the shape aesthetic is set correctly
  expect_true("shape" %in% names(plot$labels))
  expect_equal(plot$labels$shape, "type")
})

test_that("plotPCA works correctly with color and shape set to 'none'", {
  se <- create_mock_se()
  pca <- prcomp(t(assay(se)), scale = TRUE)

  plot <- plotPCA(pca, se, color = "none", shape = "none")

  # Check if the color and shape aesthetics are not set
  expect_false("colour" %in% names(plot$labels))
  expect_false("shape" %in% names(plot$labels))
})

test_that("plotPCA works correctly with color set to 'none' and shape specified", {
  se <- create_mock_se()
  pca <- prcomp(t(assay(se)), scale = TRUE)

  plot <- plotPCA(pca, se, color = "none", shape = "type")

  # Check if the shape aesthetic is set correctly
  expect_true("shape" %in% names(plot$labels))
  expect_equal(plot$labels$shape, "type")

  # Check if the color aesthetic is not set
  expect_false("colour" %in% names(plot$labels))
})

test_that("plotPCA works correctly with shape set to 'none' and color specified", {
  se <- create_mock_se()
  pca <- prcomp(t(assay(se)), scale = TRUE)

  plot <- plotPCA(pca, se, color = "group", shape = "none")

  # Check if the color aesthetic is set correctly
  expect_true("colour" %in% names(plot$labels))
  expect_equal(plot$labels$colour, "group")

  # Check if the shape aesthetic is not set
  expect_false("shape" %in% names(plot$labels))
})


######################### Tests for plotHeatmap() ##############################

test_that("plotHeatmap generates a pheatmap object", {
  se <- create_mock_se()
  
  heatmap <- plotHeatmap("Top variant", se, top = 50, title = "Top Variants Heatmap", annotationCol = NULL)
  
  # Check if the returned object is a pheatmap object
  expect_true(inherits(heatmap, "pheatmap"))
})


test_that("plotHeatmap handles 'Differentially expressed' type correctly", {
  se <- create_mock_se()
  gene_data <- as.data.frame(rowData(se))

  heatmap <- plotHeatmap("Differentially expressed", se, data = gene_data, title = "DE Genes Heatmap", annotationCol = NULL)

  # Check if the number of rows in the heatmap matches the number of genes in gene_data
  expect_equal(length(heatmap$tree_row$labels), nrow(gene_data))
})

test_that("plotHeatmap handles 'Selected time series cluster' type correctly", {
  se <- create_mock_se()
  gene_data <- rowData(se)

  heatmap <- plotHeatmap("Selected time series cluster", se, data = gene_data, title = "Clustered Genes Heatmap", annotationCol = NULL)

  # Check if the number of rows in the heatmap matches the number of unique genes in gene_data
  expect_equal(length(heatmap$tree_row$labels), length(unique(gene_data$Gene)))
})