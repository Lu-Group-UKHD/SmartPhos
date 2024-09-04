# Sample data for testing
set.seed(123)
test_matrix <- matrix(rnorm(20), nrow = 5, ncol = 4)

############################ Tests for mscale() ################################

test_that("mscale returns the same matrix when center and scale are FALSE", {
  result <- mscale(test_matrix, center = FALSE, scale = FALSE)
  expect_equal(result, test_matrix)
})

test_that("mscale correctly centers the matrix when center is TRUE and scale is FALSE", {
  result <- mscale(test_matrix, center = TRUE, scale = FALSE)
  expected <- t(apply(test_matrix, 1, function(y) y - mean(y)))
  expect_equal(result, expected)
})

test_that("mscale correctly scales the matrix when scale is TRUE and center is FALSE", {
  result <- mscale(test_matrix, center = FALSE, scale = TRUE)
  expected <- t(apply(test_matrix, 1, function(y) y / sd(y)))
  expect_equal(result, expected)
})

test_that("mscale correctly centers and scales the matrix using standard deviation", {
  result <- mscale(test_matrix, center = TRUE, scale = TRUE, useMad = FALSE)
  expected <- t(apply(test_matrix, 1, function(y) (y - mean(y)) / sd(y)))
  expect_equal(result, expected)
})

test_that("mscale correctly censors values symmetrically", {
  result <- mscale(test_matrix, center = TRUE, scale = TRUE, censor = 1)
  expected <- t(apply(test_matrix, 1, function(y) pmin(pmax((y - mean(y)) / sd(y), -1), 1)))
  expect_equal(result, expected)
})

test_that("mscale correctly censors values asymmetrically", {
  result <- mscale(test_matrix, center = TRUE, scale = TRUE, censor = c(-1, 1))
  expected <- t(apply(test_matrix, 1, function(y) pmin(pmax((y - mean(y)) / sd(y), -1), 1)))
  expect_equal(result, expected)
})

test_that("mscale works correctly when input matrix contains NA values", {
  matrix_with_na <- test_matrix
  matrix_with_na[1, 1] <- NA
  
  result <- mscale(matrix_with_na, center = TRUE, scale = TRUE, useMad = FALSE)
  expected <- t(apply(matrix_with_na, 1, function(y) (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)))
  expect_equal(result, expected)
})

test_that("mscale handles matrices with a single row correctly", {
  single_row_matrix <- matrix(rnorm(4), nrow = 1)
  result <- mscale(single_row_matrix, center = TRUE, scale = TRUE)
  expected <- (single_row_matrix - mean(single_row_matrix)) / sd(single_row_matrix)
  expect_equal(result, expected)
})


# Helper function to create a dummy SummarizedExperiment object
create_dummy_SE <- function() {
  data <- matrix(rnorm(20), nrow = 5, ncol = 4)
  colData <- DataFrame(treatment = c("TreatmentA", "TreatmentA", "Control", "Control"),
                       timepoint = c("10min", "20min", "0min", "10min"))
  colnames(data) <- paste0("Sample", 1:4)
  rownames(data) <- paste0("Gene", 1:5)
  rowData <- rownames(data)
  se <- SummarizedExperiment(assays = list(intensity = data), colData = colData, rowData = rowData)
  return(se)
}

########################### Tests for addZeroTime() ############################

test_that("addZeroTime correctly adds the zero timepoint to the specified treatment", {
  # Create a dummy SummarizedExperiment object
  data <- create_dummy_SE()
  
  # Apply addZeroTime function
  result <- addZeroTime(data, condition = "treatment", treat = "TreatmentA", zeroTreat = "Control", timeRange = c("10min", "20min"))
  
  # Expected data
  expected_assay <- assay(data)[,1:3]
  expected_colData <- DataFrame(treatment = c("TreatmentA", "TreatmentA", "TreatmentA"),
                                timepoint = c("10min", "20min", "0min"))
  rownames(expected_colData) <- paste0("Sample", 1:3)
  
  
  # Tests
  expect_equal(assay(result), expected_assay)
  expect_equal(colData(result), expected_colData)
})



# Helper function to create mock time-series data
create_mock_data <- function(rows = 10, cols = 5) {
  set.seed(123)
  mat <- matrix(rnorm(rows * cols, mean = 5, sd = 2), nrow = rows, ncol = cols)
  colnames(mat) <- paste0("Sample", 1:cols, "_", rep(1:cols), "h")
  rownames(mat) <- paste0("Gene", 1:rows)
  return(mat)
}


########################### Tests for clusterTS() ##############################

test_that("clusterTS returns a list with cluster and plot components", {
  x <- create_mock_data(rows = 10, cols = 5)
  
  result <- clusterTS(x, k = 3)
  
  expect_type(result, "list")
  expect_named(result, c("cluster", "plot"))
  expect_s3_class(result$plot, "ggplot")
  expect_s3_class(result$cluster, "tbl_df")
})

test_that("clusterTS performs clustering correctly with single condition data", {
  x <- create_mock_data(rows = 10, cols = 5)
  
  result <- clusterTS(x, k = 3)
  
  expect_true(all(result$cluster$cluster %in% paste0("cluster", 1:3)))
  expect_equal(ncol(result$cluster), 7)  # Check number of columns in cluster result
})

test_that("clusterTS performs clustering correctly with two conditions", {
  x <- create_mock_data(rows = 100, cols = 6)
  colnames(x) <- c("1h_A", "2h_A", "3h_A", "1h_B", "2h_B", "3h_B")
  
  result <- clusterTS(x, k = 3, twoCondition = TRUE)
  
  expect_true(all(result$cluster$cluster %in% paste0("cluster", 1:3)))
  expect_equal(ncol(result$cluster), 9)  # Check number of columns in cluster result
  expect_true("treatment" %in% names(result$cluster))
})

test_that("clusterTS removes rows with NA values", {
  x <- create_mock_data(rows = 10, cols = 5)
  x[1, ] <- NA
  x[2, ] <- NA
  
  result <- clusterTS(x, k = 3)
  
  expect_equal(nrow(result$cluster), 40)  # One row should be removed due to NA
})

test_that("clusterTS filters clusters based on probability cutoff", {
  x <- create_mock_data(rows = 10, cols = 5)
  
  result <- clusterTS(x, k = 3, pCut = 0.8)
  
  expect_true(all(result$cluster$prob >= 0.8))
})

test_that("clusterTS handles different time units (hours and minutes)", {
  x <- create_mock_data(rows = 10, cols = 6)
  colnames(x) <- c("1h", "30min", "2h", "45min", "3h", "10min")
  
  result <- clusterTS(x, k = 3)
  
  time_levels <- levels(result$cluster$time)
  expect_equal(time_levels, c("10min", "30min", "45min", "1h", "2h", "3h"))
})

test_that("clusterTS handles different numbers of clusters correctly", {
  x <- create_mock_data(rows = 10, cols = 5)
  
  result_k2 <- clusterTS(x, k = 2)
  result_k4 <- clusterTS(x, k = 4)
  
  expect_true(all(result_k2$cluster$cluster %in% paste0("cluster", 1:2)))
  expect_true(all(result_k4$cluster$cluster %in% paste0("cluster", 1:4)))
})

test_that("clusterTS produces reproducible results with set seed", {
  x <- create_mock_data(rows = 10, cols = 5)
  
  result1 <- clusterTS(x, k = 3)
  result2 <- clusterTS(x, k = 3)
  
  expect_equal(result1$cluster, result2$cluster)
})



######################### Tests for plotTimeSeries() ###########################

test_that("plotTimeSeries handles 'expression' type correctly", {
  # Create a mock SummarizedExperiment object
  mat <- create_mock_data(rows = 10, cols = 6)
  colData <- data.frame(timepoint=c("0h", "1h", "2h", "0h", "1h", "2h"), 
                        condition=c(rep("A", 3), rep("B", 3)), 
                        treatment=c(rep("A", 3), rep("B", 3)))
  se <- SummarizedExperiment(assays = list(intensity=mat), colData=colData)
  
  p <- plotTimeSeries(se, type = "expression", geneID = 1, symbol = "Gene 1", 
                      condition = "condition", treatment = "A", timerange = c("0h", "1h", "2h"))
  
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Gene 1")
  expect_equal(p$labels$y, "Normalized expression")
})

test_that("plotTimeSeries handles 'logFC' type correctly", {
  # Create a mock SummarizedExperiment object
  mat <- create_mock_data(rows = 10, cols = 6)
  colData <- data.frame(timepoint=c("0min", "1h", "2h", "0min", "1h", "2h"),
                        treatment=c(rep("A", 3), rep("B", 3)))
  se <- SummarizedExperiment(assays = list(intensity=mat), colData=colData)

  p <- plotTimeSeries(se, type = "logFC", geneID = "Gene1", symbol = "Gene1",
                      condition = "treatment", treatment = "A", refTreat = "B", timerange = c("0min", "1h", "2h"))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Gene1")
  expect_equal(p$labels$y, "logFC")
})

test_that("plotTimeSeries handles 'two-condition expression' type correctly", {
  # Create a mock SummarizedExperiment object
  mat <- create_mock_data(rows = 10, cols = 6)
  colData <- data.frame(timepoint=c("0min", "1h", "2h", "0min", "1h", "2h"),
                        treatment=c(rep("A", 3), rep("B", 3)))
  se <- SummarizedExperiment(assays = list(intensity=mat), colData=colData)

  p <- plotTimeSeries(se, type = "two-condition expression", geneID = "Gene1", symbol = "Gene 1",
                      condition = "treatment", treatment = "A", refTreat = "B", timerange = c("0min", "1h", "2h"))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Gene 1")
  expect_equal(p$labels$y, "Normalized expression")
})

test_that("plotTimeSeries handles 'addZero' parameter correctly", {
  # Create a mock SummarizedExperiment object
  mat <- create_mock_data(rows = 10, cols = 5)
  colData <- data.frame(timepoint=c("2h", "3h", "0min", "2h", "3h"),
                        treatment=c(rep("A", 2), rep("B", 3)))
  se <- SummarizedExperiment(assays=list(intensity=mat), colData=colData)

  p <- plotTimeSeries(se, type = "expression", geneID = "Gene1", symbol = "Gene1",
                      condition = "treatment", treatment = "A", addZero = TRUE,
                      zeroTreat = "B", timerange = c("2h", "3h"))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Gene1")
  expect_equal(p$labels$y, "Normalized expression")
})

test_that("plotTimeSeries handles subjectID correctly", {
  # Create a mock SummarizedExperiment object
  mat <- create_mock_data(rows = 10, cols = 6)
  colData <- data.frame(timepoint=c("0min", "1h", "2h", "0min", "1h", "2h"),
                        treatment=c(rep("A", 3), rep("B", 3)),
                        subjectID=rep(1:2, each=3))
  se <- SummarizedExperiment(assays=list(intensity=mat), colData=colData)

  p <- plotTimeSeries(se, type = "expression", geneID = "Gene1", symbol = "Gene 1",
                      condition = "treatment", treatment = "A", timerange = c("0min", "1h", "2h"))

  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Gene 1")
  expect_true("subjectID" %in% colnames(colData))
})