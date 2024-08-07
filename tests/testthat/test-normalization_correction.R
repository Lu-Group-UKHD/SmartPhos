########################## Tests for medianNorm() ##############################

test_that("medianNorm normalizes matrix by median", {
  x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 4, ncol = 3)
  normalized_x <- medianNorm(x, method = "median")

  expected_medians <- colMedians(x)
  adjusted_medians <- expected_medians - median(expected_medians)
  expected <- sweep(x, 2, adjusted_medians)

  expect_equal(normalized_x, expected)
})

test_that("medianNorm normalizes matrix by mean", {
  x <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12), nrow = 4, ncol = 3)
  normalized_x <- medianNorm(x, method = "mean")

  expected_means <- colMeans(x)
  adjusted_means <- expected_means - mean(expected_means)
  expected <- sweep(x, 2, adjusted_means)

  expect_equal(normalized_x, expected)
})

test_that("medianNorm handles NA values correctly with median method", {
  x <- matrix(c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10, NA, 12), nrow = 4, ncol = 3)
  normalized_x <- medianNorm(x, method = "median")

  expected_medians <- colMedians(x, na.rm = TRUE)
  adjusted_medians <- expected_medians - median(expected_medians, na.rm = TRUE)
  expected <- sweep(x, 2, adjusted_medians)

  expect_equal(normalized_x, expected)
})

test_that("medianNorm handles NA values correctly with mean method", {
  x <- matrix(c(1, 2, NA, 4, 5, 6, 7, 8, 9, 10, NA, 12), nrow = 4, ncol = 3)
  normalized_x <- medianNorm(x, method = "mean")

  expected_means <- colMeans(x, na.rm = TRUE)
  adjusted_means <- expected_means - mean(expected_means, na.rm = TRUE)
  expected <- sweep(x, 2, adjusted_means)

  expect_equal(normalized_x, expected)
})

test_that("medianNorm returns matrix with same dimensions", {
  x <- matrix(rnorm(20), nrow = 5, ncol = 4)
  normalized_x <- medianNorm(x)

  expect_equal(dim(normalized_x), dim(x))
})

test_that("medianNorm returns the same matrix when all values are the same", {
  x <- matrix(rep(5, 20), nrow = 5, ncol = 4)
  normalized_x <- medianNorm(x, method = "median")

  expected <- matrix(rep(5, 20), nrow = 5, ncol = 4)

  expect_equal(normalized_x, expected)
})



# Helper function to create a test MultiAssayExperiment object
create_test_mae <- function() {
  # Create assay data
  assay_data <- matrix(abs(rnorm(1000, 10, 2)), nrow = 100, ncol = 10)
  colnames(assay_data) <- paste0("Sample", 1:10)
  rownames(assay_data) <- paste0("Gene", 1:100)

  # Create sample metadata
  sample_data <- data.frame(sample = colnames(assay_data),
                            sampleType = rep(c("FullProteome", "Phospho"), each = 5),
                            group = rep(c("A", "B"), each = 5))
  rownames(sample_data) <- colnames(assay_data)

  # Create gene metadata
  gene_data <- data.frame(UniprotID = rownames(assay_data),
                          Gene = paste0("GeneSymbol", 1:100))
  rownames(gene_data) <- rownames(assay_data)

  # Create SummarizedExperiment object
  fpe <- SummarizedExperiment(assays = list(intensity = assay_data), colData = sample_data, rowData = gene_data)

  # Create assay data
  assay_data <- matrix(abs(rnorm(1000, 10, 5)), nrow = 100, ncol = 10)
  colnames(assay_data) <- paste0("Sample", 1:10)
  rownames(assay_data) <- paste0("Gene", 1:100)

  # Create sample metadata
  sample_data <- data.frame(sample = colnames(assay_data),
                            sampleType = rep(c("FullProteome", "Phospho"), each = 5),
                            group = rep(c("A", "B"), each = 5))
  rownames(sample_data) <- colnames(assay_data)

  # Create gene metadata
  gene_data <- data.frame(UniprotID = rownames(assay_data),
                          Gene = paste0("GeneSymbol", 1:100),
                          Residue = rep(c("S","T","Y"), c(70,15,15)),
                          Position = rep(1:100))
  rownames(gene_data) <- rownames(assay_data)

  # Create SummarizedExperiment object
  ppe <- SummarizedExperiment(assays = list(intensity = assay_data), colData = sample_data, rowData = gene_data)

  mae <- MultiAssayExperiment(list(Phosphoproteome = ppe, Proteome = fpe), colData = colData(ppe))
}


################## Tests for performCombinedNormalization() ####################

test_that("performCombinedNormalization returns matrix with correct dimensions", {
  mae <- create_test_mae()

  expect_equal(ncol(normalized_data), 5)
  expect_equal(nrow(normalized_data), 200)
})

test_that("performCombinedNormalization handles all NA rows correctly", {
  mae <- create_test_mae()

  assay(mae[["Phosphoproteome"]])[1, 1:5] <- NA

  normalized_data <- performCombinedNormalization(mae)

  expect_equal(nrow(normalized_data), 199) # One row should be removed
})


test_that("performCombinedNormalization handles different sample types correctly", {
  mae <- create_test_mae()

  colData(mae)$sampleType <- c("FullProteome", "FullProteome", "OtherType",
                               "FullProteome", "FullProteome", "Phospho",
                               "FullProteome", "Phospho", "OtherType", "FullProteome")

  normalized_data <- performCombinedNormalization(mae)

  expect_equal(ncol(normalized_data), 6)
  expect_equal(nrow(normalized_data), 200)
})



test_that("performCombinedNormalization handles non-matching sample types", {
  mae <- create_test_mae()

  colData(mae)$sampleType <- rep("OtherType", 10)

  normalized_data <- performCombinedNormalization(mae)

  expect_true(is.matrix(normalized_data))
  expect_equal(dim(normalized_data), c(0, 0))
})



