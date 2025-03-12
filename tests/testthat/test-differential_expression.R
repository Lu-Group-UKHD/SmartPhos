# Helper function to create a mock SummarizedExperiment object
create_mock_se <- function(with_subjectID = FALSE) {
  n_samples <- 10

  # Create a matrix of assay data
  assay_data <- matrix(rnorm(1000), nrow = 100, ncol = n_samples)
  rownames(assay_data) <- paste0("Gene", 1:100)
  colnames(assay_data) <- paste0("Sample", 1:10)

  # Create sample metadata
  col_data <- DataFrame(
    sample = paste0("Sample", 1:10),
    treatment = rep(c("ctr1", "EGF"), each = 5),
    timepoint = rep(c("10min", "40min"), each = 5)
  )

  # Add subjectID column if specified
  if (with_subjectID) {
    col_data$subjectID <- rep(paste0("subject", 1:(n_samples / 2)), 2)
  }

  # Create row metadata (feature data)
  gene_data <- data.frame(UniprotID = rownames(assay_data),
                          Gene = paste0("GeneSymbol", 1:100),
                          Residue = rep(c("S","T","Y"), c(70,15,15)),
                          Position = rep(1:100))
  rownames(gene_data) <- rownames(assay_data)

  # Create a SummarizedExperiment object
  se <- SummarizedExperiment(
    assays = list(Intensity = assay_data),
    colData = col_data,
    rowData = gene_data
  )

  return(se)
}


#################### Tests for performDifferentialExp() ########################

test_that("performDifferentialExp works with limma method", {
  se <- create_mock_se()

  result <- performDifferentialExp(
    se = se,
    assay = "Intensity",
    method = "limma",
    condition = "treatment",
    reference = "ctr1",
    target = "EGF"
  )

  expect_type(result, "list")
  expect_true("resDE" %in% names(result))
  expect_true("seSub" %in% names(result))
  expect_s3_class(result$resDE, "tbl_df")
  expect_s4_class(result$seSub, "SummarizedExperiment")
  expect_gt(nrow(result$resDE), 0)
})

test_that("performDifferentialExp works with ProDA method", {
  se <- create_mock_se()

  result <- performDifferentialExp(
    se = se,
    assay = "Intensity",
    method = "ProDA",
    condition = "treatment",
    reference = "ctr1",
    target = "EGF"
  )

  expect_type(result, "list")
  expect_true("resDE" %in% names(result))
  expect_true("seSub" %in% names(result))
  expect_s3_class(result$resDE, "tbl_df")
  expect_s4_class(result$seSub, "SummarizedExperiment")
  expect_gt(nrow(result$resDE), 0)
})

test_that("performDifferentialExp throws error for invalid method", {
  se <- create_mock_se()

  expect_error(
    performDifferentialExp(
      se = se,
      assay = "Intensity",
      method = "wrong_method",
      condition = "treatment",
      reference = "ctr1",
      target = "EGF"
    ),
    "Invalid method!! Provide either limma or ProDA"
  )
})

test_that("performDifferentialExp throws error for missing condition column", {
  se <- create_mock_se()

  expect_error(
    performDifferentialExp(
      se = se,
      assay = "Intensity",
      method = "limma",
      condition = "missing_column",
      reference = "ctr1",
      target = "EGF"
    ),
    "Invalid condtion! Provide condition which is part of the given Summarized Experiment object"
  )
})

test_that("performDifferentialExp correctly subsets SE object", {
  se <- create_mock_se()
  reference = "ctr1"
  target = "EGF"

  result <- performDifferentialExp(
    se = se,
    assay = "Intensity",
    method = "limma",
    condition = "treatment",
    reference = reference,
    target = target
  )

  seSub <- result$seSub
  expect_equal(length(unique(seSub$comparison)), length(reference) + length(target))  # Only reference and target groups should be present
  expect_equal(ncol(seSub), 10)  # Ensure all samples are included
})


# Helper function to create mock differential expression data
create_mock_tableDE <- function() {
  data.frame(
    ID = paste0("s", 1:10),
    log2FC = c(1.5, -2, 0.2, -0.3, 0.8, -1.7, 2.3, 0.1, -0.9, 0.5),
    pvalue = c(0.01, 0.02, 0.5, 0.6, 0.03, 0.04, 0.001, 0.8, 0.07, 0.001),
    Gene = paste0("Gene1", 1:10)
  )
}


######################### Tests for plotVolcano() ##############################

test_that("plotVolcano generates a ggplot object", {
  tableDE <- create_mock_tableDE()

  v_plot <- plotVolcano(tableDE, pFilter = 0.05, fcFilter = 0.5)

  expect_s3_class(v_plot, "ggplot")
})

test_that("plotVolcano correctly categorizes genes based on thresholds", {
  tableDE <- create_mock_tableDE()

  v_plot <- plotVolcano(tableDE, pFilter = 0.05, fcFilter = 0.5)

  dataVolcano <- ggplot2::ggplot_build(v_plot)

  up_genes <- dataVolcano$data[[8]][dataVolcano$data[[8]]$colour == "firebrick3", ]
  down_genes <- dataVolcano$data[[8]][dataVolcano$data[[8]]$colour == "navy", ]
  not_sig_genes <- dataVolcano$data[[8]][dataVolcano$data[[8]]$colour == "darkgrey", ]

  expect_equal(nrow(up_genes), 4)  # Expected number of upregulated genes
  expect_equal(nrow(down_genes), 2)  # Expected number of downregulated genes
  expect_equal(nrow(not_sig_genes), 4)  # Expected number of non-significant genes
})

test_that("plotVolcano handles non-default column names correctly", {
  tableDE <- create_mock_tableDE()

  names(tableDE) <- c("ID", "log2FoldChange", "pvalue", "Gene")
  expect_error(plotVolcano(tableDE), "column 'log2FC' not found")

  names(tableDE) <- c("ids", "log2FC", "pvalue", "Gene")
  expect_error(plotVolcano(tableDE), "column 'ID' not found")

  names(tableDE) <- c("ID", "log2FC", "p_value", "Gene")
  expect_error(plotVolcano(tableDE), "column 'pvalue' not found")

  names(tableDE) <- c("ID", "log2FC", "pvalue", "GeneName")
  expect_error(plotVolcano(tableDE), "column 'Gene' not found")

})



########################### Tests for intensityBoxPlot() ################################

test_that("intensityBoxPlot generates a ggplot object without subjectID", {
  se <- create_mock_se(with_subjectID = FALSE)
  result <- performDifferentialExp(
    se = se,
    assay = "Intensity",
    method = "limma",
    condition = "treatment",
    reference = "ctr1",
    target = "EGF"
  )
  p <- intensityBoxPlot(result$seSub, id = "Gene1", symbol = "GeneSymbol1")

  expect_s3_class(p, "ggplot")
})

test_that("intensityBoxPlot generates a ggplot object with subjectID", {
  se <- create_mock_se(with_subjectID = TRUE)
  result <- performDifferentialExp(
    se = se,
    assay = "Intensity",
    method = "limma",
    condition = "treatment",
    reference = "ctr1",
    target = "EGF"
  )
  p <- intensityBoxPlot(result$seSub, id = "Gene1", symbol = "GeneSymbol1")

  expect_s3_class(p, "ggplot")
})

test_that("intensityBoxPlot handles missing gene/feature ID gracefully", {
  se <- create_mock_se(with_subjectID = FALSE)
  result <- performDifferentialExp(
    se = se,
    assay = "Intensity",
    method = "limma",
    condition = "treatment",
    reference = "ctr1",
    target = "EGF"
  )

  expect_error(intensityBoxPlot(result$seSub, id = "nonexistent_gene", symbol = "Nonexistent Gene"),
               "subscript out of bounds")
})

test_that("intensityBoxPlot includes correct elements in the plot with subjectID", {
  se <- create_mock_se(with_subjectID = TRUE)
  result <- performDifferentialExp(
    se = se,
    assay = "Intensity",
    method = "limma",
    condition = "treatment",
    reference = "ctr1",
    target = "EGF"
  )
  p <- intensityBoxPlot(result$seSub, id = "Gene1", symbol = "GeneSymbol1")
  p_data <- ggplot2::ggplot_build(p)$data

  # Expect a boxplot layer, points layer, and line layer
  expect_equal(length(p_data), 3)

  # Check that the boxplot, points, and lines are based on the correct groups and subjects
  expect_equal(unique(p_data[[1]]$x), c(1, 2))  # Group positions on x-axis
  expect_equal(length(unique(p_data[[2]]$group)), 2)  # For 2 groups in comparison
})

test_that("intensityBoxPlot handles different comparison groups", {
  se <- create_mock_se(with_subjectID = FALSE)

  # Modify comparison groups
  se$comparison <- rep(c("A", "B", "C"), length.out = ncol(se))

  p <- intensityBoxPlot(se, id = "Gene1", symbol = "GeneSymbol1")
  p_data <- ggplot2::ggplot_build(p)$data

  # Check that there are 3 groups plotted
  expect_equal(length(unique(p_data[[1]]$x)), 3)
})
