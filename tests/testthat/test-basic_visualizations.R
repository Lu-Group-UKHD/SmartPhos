######################## Tests for plotMissing() ###############################

# Helper function to create a mock SummarizedExperiment
create_mock_se <- function() {
  # Create sample data
  assay_data <- matrix(c(1, NA, 3, 4, 5, NA, 7, 8, 9, 10, NA, 12, 13, 14, 15), nrow = 3, ncol = 5)
  colnames(assay_data) <- c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5")
  
  # Create sample annotations
  sample_data <- data.frame(sample = colnames(assay_data))
  rownames(sample_data) <- colnames(assay_data)
  
  # Create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(counts = assay_data), colData = sample_data)
  
  return(se)
}

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