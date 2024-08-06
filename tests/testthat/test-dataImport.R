##################### Tests for readExperimentDIA() ###########################

# Mock data for testing
file1 <- system.file("extdata", "phosDIA_1.xls", package = "SmartPhos")
file2 <- system.file("extdata", "proteomeDIA_1.xls", package = "SmartPhos")

fileTable <- data.frame(
  type = c("phosphoproteome", "proteome", "proteome"),
  fileName = c(file1, file2, file2),
  id = c("Sample_1", "sample1", "sample2"),
  outputID = c("s1", "s2", "s3")
)


# Unit tests
test_that("readExperimentDIA processes both phosphoproteome and proteome data correctly", {
  result <- readExperimentDIA(fileTable, localProbCut = 0.75, annotation_col = c("id"), onlyReviewed = FALSE, normalizeByProtein = FALSE)

  expect_s4_class(result, "MultiAssayExperiment")
  expect_true("Phosphoproteome" %in% names(experiments(result)))
  expect_true("Proteome" %in% names(experiments(result)))
  expect_equal(nrow(colData(result)), 3)
  expect_equal(ncol(experiments(result)$Phosphoproteome), 1)
  expect_equal(ncol(experiments(result)$Proteome), 2)
  expect_equal(rownames(colData(result)), c("s1", "s2", "s3"))
  expect_equal(colData(result)$sample, c("s1", "s2", "s3"))
})

test_that("readExperimentDIA handles missing phosphoproteome data", {

  fileTable <- data.frame(
    type = c("proteome", "proteome"),
    fileName = c(file2, file2),
    id = c("sample1", "sample2")
  )

  result <- readExperimentDIA(fileTable, localProbCut = 0.75, annotation_col = c("id"), normalizeByProtein = FALSE)

  expect_s4_class(result, "MultiAssayExperiment")
  expect_false("Phosphoproteome" %in% names(experiments(result)))
  expect_true("Proteome" %in% names(experiments(result)))
  expect_equal(nrow(colData(result)), 2)
  expect_equal(ncol(experiments(result)$Proteome), 2)
})

test_that("readExperimentDIA handles missing proteome data", {

  fileTable <- data.frame(
    type = c("phosphoproteome"),
    fileName = c(file1),
    id = c("Sample_1")
  )

  result <- readExperimentDIA(fileTable, localProbCut = 0.75, annotation_col = c("id"), onlyReviewed = FALSE, normalizeByProtein = FALSE)

  expect_s4_class(result, "MultiAssayExperiment")
  expect_true("Phosphoproteome" %in% names(experiments(result)))
  expect_false("Proteome" %in% names(experiments(result)))
  expect_equal(nrow(colData(result)), 1)
  expect_equal(ncol(experiments(result)$Phosphoproteome), 1)
})


####################### Tests for readExperiment() #############################


# Mock data for testing
file1 <- system.file("extdata", "phosDDA_1.xls", package = "SmartPhos")
file2 <- system.file("extdata", "proteomeDDA_1.xls", package = "SmartPhos")

fileTable <- data.frame(
  type = c("phosphoproteome", "proteome"),
  fileName = c(file1, file2),
  sample = c("Sample1", "sample1"),
  id = c("s1", "s2")
)

# Unit tests
test_that("readExperiment processes both phosphoproteomic and proteomic data correctly", {

  result <- readExperiment(fileTable, localProbCut = 0.75, scoreDiffCut = 5, fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE, annotation_col = c("id"))

  expect_s4_class(result, "MultiAssayExperiment")
  expect_true("Phosphoproteome" %in% names(experiments(result)))
  expect_true("Proteome" %in% names(experiments(result)))
  expect_equal(nrow(colData(result)), 2)
  expect_equal(ncol(experiments(result)$Phosphoproteome), 1)
  expect_equal(ncol(experiments(result)$Proteome), 1)
})

test_that("readExperiment handles missing phosphoproteomic data", {
  fileTable <- data.frame(
    type = c("proteome"),
    fileName = c(file2),
    sample = c("sample1"),
    id = c("s1")
  )

  result <- readExperiment(fileTable, localProbCut = 0.75, scoreDiffCut = 5, fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE, annotation_col = c())

  expect_s4_class(result, "MultiAssayExperiment")
  expect_false("Phosphoproteome" %in% names(experiments(result)))
  expect_true("Proteome" %in% names(experiments(result)))
  expect_equal(nrow(colData(result)), 1)
  expect_equal(ncol(experiments(result)$Proteome), 1)
})

test_that("readExperiment handles missing proteomic data", {
  fileTable <- data.frame(
    type = c("phosphoproteome"),
    fileName = c(file1),
    sample = c("Sample1"),
    id = c("s1")
  )

  result <- readExperiment(fileTable, localProbCut = 0.75, scoreDiffCut = 5, fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE, annotation_col = c())

  expect_s4_class(result, "MultiAssayExperiment")
  expect_true("Phosphoproteome" %in% names(experiments(result)))
  expect_false("Proteome" %in% names(experiments(result)))
  expect_equal(nrow(colData(result)), 1)
  expect_equal(ncol(experiments(result)$Phosphoproteome), 1)
})

test_that("readExperiment constructs sample annotations correctly", {
  fileTable2 <- data.frame(
    type = c("phosphoproteome", "proteome"),
    fileName = c(file1, file2),
    sample = c("Sample1", "sample1"),
    id = c("s1", "s2"),
    batch = c("batch1", "batch2"),
    condition = c("control", "treated")
  )

  result <- readExperiment(fileTable2, localProbCut = 0.75, scoreDiffCut = 5, fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE, annotation_col = c("batch", "condition"))

  expect_equal(rownames(colData(result)), c("s1", "s2"))
  expect_equal(colData(result)$sample, c("Sample1", "sample1"))
  expect_equal(colData(result)$batch, c("batch1", "batch2"))
  expect_equal(colData(result)$condition, c("control", "treated"))
})


##################### Tests for normByFullProteome() ###########################

# Helper function to create a mock MultiAssayExperiment
create_mock_mae <- function() {
  # Create sample data
  phospho_data <- matrix(runif(20), nrow=4, ncol=5)
  proteo_data <- matrix(runif(20, 1, 2), nrow=4, ncol=5)
  
  # Set column names to match sample names
  colnames(phospho_data) <- colnames(proteo_data) <- c("S1", "S2", "S3", "S4", "S5")
  
  # Create sample annotations
  sample_data <- data.frame(sampleName = c("S1", "S2", "S3", "S4", "S5"),
                            sampleType = c("FullProteome", "FullProteome", "Phosphoproteome", "Phosphoproteome", "Phosphoproteome"))
  rownames(sample_data) <- colnames(phospho_data)
  
  # Create row annotations
  row_data <- data.frame(UniprotID = paste0("P", 1:4))
  
  # Create SummarizedExperiment objects
  ppe <- SummarizedExperiment(assays = list(counts = phospho_data), colData = sample_data, rowData = row_data)
  fpe <- SummarizedExperiment(assays = list(counts = proteo_data), colData = sample_data, rowData = row_data)
  
  # Create MultiAssayExperiment
  mae <- MultiAssayExperiment(experiments = list(Phosphoproteome = ppe, Proteome = fpe), colData = sample_data)
  
  return(mae)
}

test_that("normByFullProteome checks for required assays", {
  mae <- create_mock_mae()

  # Remove Phosphoproteome assay
  mae_no_phospho <- mae
  mae_no_phospho[ , , "Phosphoproteome"] <- NULL

  expect_error(normByFullProteome(mae_no_phospho), "Both Phosphoproteome and Proteome assays should be present in the MultiAssayExperiment object")

  # Remove Proteome assay
  mae_no_proteome <- mae
  mae_no_proteome[ , , "Proteome"] <- NULL

  expect_error(normByFullProteome(mae_no_proteome), "Both Phosphoproteome and Proteome assays should be present in the MultiAssayExperiment object")
})

test_that("normByFullProteome handles missing samples correctly", {
  mae <- create_mock_mae()
  colData(mae)$sampleType <- "Phosphoproteome"
  
  expect_error(normByFullProteome(mae), "Proteome assay for the unenriched samples i.e., sampleType with FullProteome should be present")
})

