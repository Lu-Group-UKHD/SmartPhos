# Helper function to create a mock SummarizedExperiment object
create_mock_data <- function() {
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
                          Gene = paste0("GeneSymbol", 1:100), 
                          Residue = rep(c("S","T","Y"), c(70,15,15)),
                          Position = rep(1:100))
  rownames(gene_data) <- rownames(assay_data)
  
  # Create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(intensity = assay_data), colData = sample_data, rowData = gene_data)
  
  return(se)
}


##################### Tests for preprocessProteome() ###########################

test_that("preprocessProteome returns a SummarizedExperiment object", {
  seData <- create_mock_data()

  processedData <- preprocessProteome(seData)

  # Check if the returned object is a SummarizedExperiment
  expect_true(inherits(processedData, "SummarizedExperiment"))
})

test_that("preprocessProteome filters based on sample type", {
  seData <- create_mock_data()

  processedData <- preprocessProteome(seData, getPP = TRUE)

  # Check if only PP samples are included
  expect_true(all(colData(processedData)$sampleType %in% c("Phospho", "PP")))
})

test_that("preprocessProteome removes specified outliers", {
  seData <- create_mock_data()

  processedData <- preprocessProteome(seData, removeOutlier = "Sample1")

  # Check if Sample1 is removed
  expect_false("Sample1" %in% colnames(processedData))
})

test_that("preprocessProteome applies filters correctly", {
  seData <- create_mock_data()

  processedData <- preprocessProteome(seData, filterList = list(group = "A"))

  # Check if only samples in group A are included
  expect_true(all(colData(processedData)$group == "A"))
})

test_that("preprocessProteome performs log2 transformation", {
  seData <- create_mock_data()

  processedData <- preprocessProteome(seData, transform = "log2")

  # Check if the data is log2 transformed
  expect_true(all(assay(processedData) <= log2(max(assay(seData))) &
                    assay(processedData) >= log2(min(assay(seData)))))
})

test_that("preprocessProteome performs imputation", {
  seData <- create_mock_data()

  # Introduce some missing values
  assay(seData)[1:10, 1:2] <- NA

  processedData <- preprocessProteome(seData, impute = "QRILC")

  # Check if missing values are imputed
  expect_true(all(!is.na(assays(processedData)[["imputed"]][1:10, 1:2])))
})



######################## Tests for preprocessPhos() ############################

test_that("preprocessPhos returns a SummarizedExperiment object", {
  seData <- create_mock_data()
  
  processedData <- preprocessPhos(seData)
  
  # Check if the returned object is a SummarizedExperiment
  expect_true(inherits(processedData, "SummarizedExperiment"))
})

test_that("preprocessPhos filters based on sample type", {
  seData <- create_mock_data()

  processedData <- preprocessPhos(seData, getFP = TRUE)

  # Check if only FP samples are included
  expect_true(all(colData(processedData)$sampleType %in% c("FullProteome", "FP")))
})

test_that("preprocessPhos removes specified outliers", {
  seData <- create_mock_data()

  processedData <- preprocessPhos(seData, removeOutlier = "Sample1")

  # Check if Sample1 is removed
  expect_false("Sample1" %in% colnames(processedData))
})

test_that("preprocessPhos applies filters correctly", {
  seData <- create_mock_data()

  processedData <- preprocessPhos(seData, filterList = list(group = "B"))

  # Check if only samples in group A are included
  expect_true(all(colData(processedData)$group == "B"))
})

test_that("preprocessPhos performs log2 transformation", {
  seData <- create_mock_data()

  processedData <- preprocessPhos(seData, transform = "log2")

  # Check if the data is log2 transformed
  expect_true(all(assay(processedData) <= log2(max(assay(seData))) &
                    assay(processedData) >= log2(min(assay(seData)))))
})

test_that("preprocessPhos performs imputation", {
  seData <- create_mock_data()

  # Introduce some missing values
  assay(seData)[1:10, 1:2] <- NA

  processedData <- preprocessPhos(seData, impute = "QRILC")

  # Check if missing values are imputed
  expect_true(all(!is.na(assays(processedData)[["imputed"]][1:10, 1:2])))
})