# # Mock data for testing
# file1 <- system.file("extdata", "phosDIA_1.xls", package = "SmartPhos")
# file2 <- system.file("extdata", "proteomeDIA_1.xls", package = "SmartPhos")
# 
# fileTable <- data.frame(
#   type = c("phosphoproteome", "proteome", "proteome"),
#   fileName = c(file1, file2, file2),
#   id = c("Sample_1", "sample1", "sample2"),
#   outputID = c("s1", "s2", "s3")
# )
# 
# 
# # Unit tests
# test_that("readExperimentDIA processes both phosphoproteome and proteome data correctly", {
#   result <- readExperimentDIA(fileTable, localProbCut = 0.75, annotation_col = c("id"), onlyReviewed = FALSE, normalizeByProtein = FALSE)
#   
#   expect_s4_class(result, "MultiAssayExperiment")
#   expect_true("Phosphoproteome" %in% names(experiments(result)))
#   expect_true("Proteome" %in% names(experiments(result)))
#   expect_equal(nrow(colData(result)), 3)
#   expect_equal(ncol(experiments(result)$Phosphoproteome), 1)
#   expect_equal(ncol(experiments(result)$Proteome), 2)
#   expect_equal(rownames(colData(result)), c("s1", "s2", "s3"))
#   expect_equal(colData(result)$sample, c("s1", "s2", "s3"))
# })
# 
# test_that("readExperimentDIA handles missing phosphoproteome data", {
#   
#   fileTable <- data.frame(
#     type = c("proteome", "proteome"),
#     fileName = c(file2, file2),
#     id = c("sample1", "sample2")
#   )
# 
#   result <- readExperimentDIA(fileTable, localProbCut = 0.75, annotation_col = c("id"), normalizeByProtein = FALSE)
# 
#   expect_s4_class(result, "MultiAssayExperiment")
#   expect_false("Phosphoproteome" %in% names(experiments(result)))
#   expect_true("Proteome" %in% names(experiments(result)))
#   expect_equal(nrow(colData(result)), 2)
#   expect_equal(ncol(experiments(result)$Proteome), 2)
# })
# 
# test_that("readExperimentDIA handles missing proteome data", {
#   
#   fileTable <- data.frame(
#     type = c("phosphoproteome"),
#     fileName = c(file1),
#     id = c("Sample_1")
#   )
# 
#   result <- readExperimentDIA(fileTable, localProbCut = 0.75, annotation_col = c("id"), onlyReviewed = FALSE, normalizeByProtein = FALSE)
# 
#   expect_s4_class(result, "MultiAssayExperiment")
#   expect_true("Phosphoproteome" %in% names(experiments(result)))
#   expect_false("Proteome" %in% names(experiments(result)))
#   expect_equal(nrow(colData(result)), 1)
#   expect_equal(ncol(experiments(result)$Phosphoproteome), 1)
# })




# Mock data for testing
file1 <- system.file("extdata", "phosDDA_1.xls", package = "SmartPhos")
file2 <- system.file("extdata", "proteomeDDA_1.xls", package = "SmartPhos")

fileTable <- data.frame(
  type = c("phosphoproteome", "proteome"),
  fileName = c(file1, file2),
  id = c("Sample1", "sample1"),
  outputID = c("s1", "s2")
)

# Unit tests
test_that("readExperiment processes both phosphoproteomic and proteomic data correctly", {
  
  print(fileTable)
  result <- readExperiment(fileTable, localProbCut = 0.75, scoreDiffCut = 5, fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE, annotation_col = c("id"))
  
  expect_s4_class(result, "MultiAssayExperiment")
  # expect_true("Phosphoproteome" %in% names(experiments(result)))
  # expect_true("Proteome" %in% names(experiments(result)))
  # expect_equal(nrow(colData(result)), 2)
  # expect_equal(ncol(experiments(result)$Phosphoproteome), 2)
  # expect_equal(ncol(experiments(result)$Proteome), 2)
})

# test_that("readExperiment handles missing phosphoproteomic data", {
#   readPhosphoExperiment <- function(...) NULL
#   
#   result <- readExperiment(fileTable, localProbCut = 0.75, scoreDiffCut = 5, fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE, annotation_col = c())
#   
#   expect_s4_class(result, "MultiAssayExperiment")
#   expect_false("Phosphoproteome" %in% names(experiments(result)))
#   expect_true("Proteome" %in% names(experiments(result)))
#   expect_equal(nrow(colData(result)), 2)
#   expect_equal(ncol(experiments(result)$Proteome), 2)
# })
# 
# test_that("readExperiment handles missing proteomic data", {
#   readProteomeExperiment <- function(...) NULL
#   
#   result <- readExperiment(fileTable, localProbCut = 0.75, scoreDiffCut = 5, fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE, annotation_col = c())
#   
#   expect_s4_class(result, "MultiAssayExperiment")
#   expect_true("Phosphoproteome" %in% names(experiments(result)))
#   expect_false("Proteome" %in% names(experiments(result)))
#   expect_equal(nrow(colData(result)), 2)
#   expect_equal(ncol(experiments(result)$Phosphoproteome), 2)
# })
# 
# test_that("readExperiment constructs sample annotations correctly", {
#   fileTable2 <- data.frame(
#     id = c("sample1", "sample2"),
#     sample = c("A", "B"),
#     type = c("phosphoproteome", "proteome"),
#     batch = c("batch1", "batch2")
#   )
#   
#   result <- readExperiment(fileTable2, localProbCut = 0.75, scoreDiffCut = 5, fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE, annotation_col = c("batch"))
#   
#   expect_equal(rownames(colData(result)), c("sample1", "sample2"))
#   expect_equal(colData(result)$sample, c("A", "B"))
#   expect_equal(colData(result)$batch, c("batch1", "batch2"))
# })
# 
# test_that("readExperiment uses additional annotation columns correctly", {
#   fileTable3 <- data.frame(
#     id = c("sample1", "sample2"),
#     sample = c("A", "B"),
#     type = c("phosphoproteome", "proteome"),
#     condition = c("control", "treated")
#   )
#   
#   result <- readExperiment(fileTable3, localProbCut = 0.75, scoreDiffCut = 5, fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE, annotation_col = c("condition"))
#   
#   expect_equal(rownames(colData(result)), c("sample1", "sample2"))
#   expect_equal(colData(result)$sample, c("A", "B"))
#   expect_equal(colData(result)$condition, c("control", "treated"))
# })