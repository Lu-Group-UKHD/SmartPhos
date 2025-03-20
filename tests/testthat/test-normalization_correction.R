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
create_test_mae <- function(adjusted = FALSE) {
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
  fpe <- SummarizedExperiment(assays = list(Intensity = assay_data),
                              colData = sample_data, rowData = gene_data)

  # Create assay data
  assay_data <- matrix(abs(rnorm(1000, 10, 5)), nrow = 100, ncol = 10)
  colnames(assay_data) <- paste0("Sample", 1:10)
  rownames(assay_data) <- paste0("Gene", 1:100)

  # Create sample metadata
  sample_data <- data.frame(sample = colnames(assay_data),
                            sampleType = rep(c("FullProteome", "Phospho"), each = 5),
                            sampleName = c("s1","s2","s3","s4","s5","s1",
                                           "s2","s3","s4","s5"),
                            group = rep(c("A", "B"), each = 5))
  rownames(sample_data) <- colnames(assay_data)

  # Create gene metadata
  gene_data <- data.frame(UniprotID = rownames(assay_data),
                          Gene = paste0("GeneSymbol", 1:100),
                          Residue = rep(c("S","T","Y"), c(70,15,15)),
                          Position = rep(1:100))
  rownames(gene_data) <- rownames(assay_data)

  # Create SummarizedExperiment object
  if (adjusted) {
    sample_data$adjustFactorPP <- runif(10, 0.8, 1.2)
    assay_data_adjusted = assay_data * runif(10, 0.8, 1.2)
    ppe <- SummarizedExperiment(assays = list(Intensity = assay_data,
                                              Intensity_adjusted = assay_data_adjusted),
                                colData = sample_data, rowData = gene_data)

  } else {
    ppe <- SummarizedExperiment(assays = list(Intensity = assay_data),
                                colData = sample_data, rowData = gene_data)
  }

  mae <- MultiAssayExperiment(list(Phosphoproteome = ppe, Proteome = fpe),
                              colData = colData(ppe))
}


################## Tests for performCombinedNormalization() ####################

test_that("performCombinedNormalization returns matrix with correct dimensions", {
  mae <- create_test_mae()

  normalized_data <- performCombinedNormalization(mae)

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
                               "FullProteome", "Phospho", "OtherType",
                               "FullProteome")

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


######################## Tests for getRatioMatrix() ############################

test_that("getRatioMatrix returns matrix with correct dimensions", {
  mae <- create_test_mae()
  ratio_matrix <- getRatioMatrix(mae)

  expect_true(is.matrix(ratio_matrix))
  expect_equal(ncol(ratio_matrix), 5)
  expect_equal(nrow(ratio_matrix), 100)
})

test_that("getRatioMatrix handles all NA rows correctly", {
  mae <- create_test_mae()

  assay(mae[["Phosphoproteome"]])[1, 1:5] <- NA

  ratio_matrix <- getRatioMatrix(mae)

  expect_equal(nrow(ratio_matrix), 99)
})


test_that("getRatioMatrix handles adjusted phosphoproteome data correctly", {
  mae <- create_test_mae()

  # Adding an adjusted phosphoproteome assay
  assays(mae[["Phosphoproteome"]])$Intensity_adjusted <- assay(mae[["Phosphoproteome"]]) * 1.2
  ratio_matrix <- getRatioMatrix(mae, getAdjustedPP = TRUE)

  expect_equal(ncol(ratio_matrix), 5)
  expect_equal(nrow(ratio_matrix), 100)
})


test_that("getRatioMatrix handles non-matching sample names", {
  mae <- create_test_mae()

  colData(mae)$sampleName <- c(rep("OtherType1", 5), rep("OtherType2", 5))

  ratio_matrix <- getRatioMatrix(mae)

  expect_true(is.matrix(ratio_matrix))
  expect_equal(dim(ratio_matrix), c(0, 0))
})


######################### Tests for plotLogRatio() #############################

# Test if plotLogRatio returns a ggplot object
test_that("plotLogRatio returns a ggplot object", {
  mae <- create_test_mae()
  plot <- plotLogRatio(mae, normalization = FALSE)

  expect_s3_class(plot, "ggplot")
})

# Test if plotLogRatio works without normalization
test_that("plotLogRatio works without normalization", {
  mae <- create_test_mae()
  plot <- plotLogRatio(mae, normalization = FALSE)

  # Check if plot has the correct title
  expect_equal(plot$labels$title,
               "Boxplot of Phospho/FullProteome Ratio(un-adjusted)")

  # Check if the x and y labels are correct
  expect_equal(plot$labels$x, "sample")
  expect_equal(plot$labels$y, "log2(ratio)")

  # Check if the plot has a horizontal line (geom_hline)
  expect_true(any(sapply(plot$layers, function(x) inherits(x$geom, "GeomHline"))))
})

# Test if plotLogRatio works with normalization
test_that("plotLogRatio works with normalization", {
  mae <- create_test_mae()

  # Add a dummy performCombinedNormalization function
  performCombinedNormalization <- function(maeData) {
    combined_data <- rbind(assay(maeData[["Proteome"]]),
                           assay(maeData[["Phosphoproteome"]]))
    return(combined_data)
  }

  plot <- plotLogRatio(mae, normalization = TRUE)

  expect_s3_class(plot, "ggplot")

  # Check if plot has the correct title
  expect_equal(plot$labels$title,
               "Boxplot of Phospho/FullProteome Ratio(un-adjusted)")
})

# Test if plotLogRatio handles missing data correctly
test_that("plotLogRatio handles missing data correctly", {
  mae <- create_test_mae()

  # Introduce some NA values
  assay(mae[["Proteome"]])[1, 1] <- NA
  assay(mae[["Phosphoproteome"]])[2, 2] <- NA

  plot <- plotLogRatio(mae, normalization = FALSE)

  # Check if plot is created successfully
  expect_s3_class(plot, "ggplot")
})



######################### Tests for checkRatioMat() ############################


# Test if checkRatioMat returns an empty vector when all samples meet the criteria
test_that("checkRatioMat returns NULL when all samples meet the criteria", {
  ratioMat <- matrix(c(1, 1, 1, NA, 1, 1, 1, 1, 1, 1, 1, 1), nrow = 4, ncol = 3)
  colnames(ratioMat) <- c("Sample1", "Sample2", "Sample3")

  excluded_samples <- checkRatioMat(ratioMat, minOverlap = 2)

  expect_equal(excluded_samples, NULL)
})


# Test if checkRatioMat handles empty ratioMat correctly
test_that("checkRatioMat handles empty ratioMat correctly", {
  ratioMat <- matrix(numeric(0), nrow = 0, ncol = 0)

  excluded_samples <- checkRatioMat(ratioMat, minOverlap = 3)

  expect_equal(excluded_samples, NULL)
})



####################### Tests for runPhosphoAdjustment() #######################

# Test if runPhosphoAdjustment returns a MultiAssayExperiment object
test_that("runPhosphoAdjustment returns a MultiAssayExperiment object", {
  mae <- create_test_mae()

  adjusted_mae <- runPhosphoAdjustment(mae, normalization = FALSE,
                                       minOverlap = 3, completeness = 0.5,
                                       ncore = 1)

  expect_s4_class(adjusted_mae, "MultiAssayExperiment")
})

# Test if runPhosphoAdjustment correctly adds adjusted intensity assay
test_that("runPhosphoAdjustment adds adjusted intensity assay", {
  mae <- create_test_mae()

  adjusted_mae <- runPhosphoAdjustment(mae, normalization = FALSE,
                                       minOverlap = 3, completeness = 0.5,
                                       ncore = 1)

  expect_true("Intensity_adjusted" %in% assayNames(experiments(adjusted_mae)[["Phosphoproteome"]]))
})

# Test if runPhosphoAdjustment correctly adjusts intensities
test_that("runPhosphoAdjustment correctly adjusts phosphoproteome intensities", {
  mae <- create_test_mae()

  # Perform the adjustment
  adjusted_mae <- runPhosphoAdjustment(mae, normalization = FALSE,
                                       minOverlap = 3, completeness = 0.5,
                                       ncore = 1)

  # Extract original and adjusted intensities
  original_intensities <- assay(mae[["Phosphoproteome"]])
  adjusted_intensities <- assays(adjusted_mae[["Phosphoproteome"]])[["Intensity_adjusted"]]

  # Expect that adjusted intensities differ from original, indicating adjustment has occurred
  expect_false(identical(original_intensities, adjusted_intensities))
})

# Test if runPhosphoAdjustment handles different values of minOverlap correctly
test_that("runPhosphoAdjustment handles minOverlap parameter correctly", {
  mae <- create_test_mae()

  # Adjust with minOverlap = 3
  adjusted_mae_min3 <- runPhosphoAdjustment(mae, normalization = FALSE,
                                            minOverlap = 3, completeness = 0.5,
                                            ncore = 1)
  samples_min3 <- colnames(assays(adjusted_mae_min3[["Phosphoproteome"]])[["Intensity_adjusted"]])

  # Adjust with minOverlap = 5
  adjusted_mae_min5 <- runPhosphoAdjustment(mae, normalization = FALSE,
                                            minOverlap = 5, completeness = 0.5,
                                            ncore = 1)
  samples_min5 <- colnames(assays(adjusted_mae_min5[["Phosphoproteome"]])[["Intensity_adjusted"]])

  # Expect that the number of samples may differ due to stricter overlap requirements
  expect_true(length(samples_min5) <= length(samples_min3))
})

# Test if runPhosphoAdjustment handles different values of completeness correctly
test_that("runPhosphoAdjustment handles completeness parameter correctly", {
  mae <- create_test_mae()

  # Adjust with completeness = 0.5
  adjusted_mae_comp05 <- runPhosphoAdjustment(mae, normalization = FALSE,
                                              minOverlap = 3,
                                              completeness = 0.5, ncore = 1)
  features_comp05 <- rownames(assays(adjusted_mae_comp05[["Phosphoproteome"]])[["Intensity_adjusted"]])

  # Adjust with completeness = 1.0
  adjusted_mae_comp10 <- runPhosphoAdjustment(mae, normalization = FALSE,
                                              minOverlap = 3,
                                              completeness = 1.0, ncore = 1)
  features_comp10 <- rownames(assays(adjusted_mae_comp10[["Phosphoproteome"]])[["Intensity_adjusted"]])

  # Expect that the number of features may differ due to stricter completeness requirements
  expect_true(length(features_comp10) <= length(features_comp05))
})

# Test if runPhosphoAdjustment correctly handles ncore > 1 (parallel processing)
test_that("runPhosphoAdjustment handles parallel processing with ncore > 1", {
  mae <- create_test_mae()

  # Adjust with parallel processing
  adjusted_mae_parallel <- runPhosphoAdjustment(mae, normalization = FALSE,
                                                minOverlap = 3,
                                                completeness = 0.5, ncore = 2)

  expect_s4_class(adjusted_mae_parallel, "MultiAssayExperiment")
})



####################### Tests for plotAdjustmentResults() ######################

# Test if plotAdjustmentResults returns a list
test_that("plotAdjustmentResults returns a list", {
  mae <- create_test_mae(adjusted = TRUE)

  plots <- plotAdjustmentResults(mae, normalization = FALSE)

  expect_type(plots, "list")
})

# Test if plotAdjustmentResults returns all required plots when adjustment=TRUE
test_that("plotAdjustmentResults returns all required plots", {
  mae <- create_test_mae(adjusted = TRUE)

  plots <- plotAdjustmentResults(mae, normalization = FALSE)

  expect_true("ratioTrendPlot" %in% names(plots))
  expect_true("ratioBoxplot" %in% names(plots))
  expect_true("ppBoxplot" %in% names(plots))

  expect_s3_class(plots$ratioTrendPlot, "gg")
  expect_s3_class(plots$ratioBoxplot, "gg")
  expect_s3_class(plots$ppBoxplot, "gg")
})

# Test if plotAdjustmentResults throws an error when adjustment factor is not present
test_that("plotAdjustmentResults throws error when adjustment factor is missing", {
  mae <- create_test_mae(adjusted = FALSE)

  expect_error(plotAdjustmentResults(mae, normalization = FALSE),
               "Phosphorylation measurments have not been adjusted yet. Please perform normalization adjustment using calcAdjustFacotr function first")
})

# Test if plotAdjustmentResults generates ratioTrendPlot correctly when all
# features are present in all samples
test_that("plotAdjustmentResults generates ratioTrendPlot correctly", {
  mae <- create_test_mae(adjusted = TRUE)

  plots <- plotAdjustmentResults(mae, normalization = FALSE)

  expect_s3_class(plots$ratioTrendPlot, "gg")
})

# Test if plotAdjustmentResults works correctly with normalization set to TRUE
test_that("plotAdjustmentResults works with normalization = TRUE", {
  mae <- create_test_mae(adjusted = TRUE)

  plots <- plotAdjustmentResults(mae, normalization = TRUE)

  expect_s3_class(plots$ratioTrendPlot, "gg")
  expect_s3_class(plots$ratioBoxplot, "gg")
  expect_s3_class(plots$ppBoxplot, "gg")
})

# Test if plotAdjustmentResults plots are not empty
test_that("plotAdjustmentResults plots are not empty", {
  mae <- create_test_mae(adjusted = TRUE)

  plots <- plotAdjustmentResults(mae, normalization = FALSE)

  # Check that plots are not empty
  expect_gt(length(ggplot2::ggplot_build(plots$ratioTrendPlot)$data), 0)
  expect_gt(length(ggplot2::ggplot_build(plots$ratioBoxplot)$data), 0)
  expect_gt(length(ggplot2::ggplot_build(plots$ppBoxplot)$data), 0)
})
