

test_that("readOneProteomDIA returns a data.table", {
  inputTab <- data.table(
    Sample1.PG.Quantity. = c("1000", "Filtered", "2000"),
    PG.ProteinGroups = c("P1", "P2", "P3"),
    PG.Genes = c("Gene1", "Gene2", "Gene3")
  )
  
  result <- readOneProteomDIA(inputTab, "Sample1")
  
  expect_s3_class(result, "data.table")
})

test_that("readOneProteomDIA correctly processes and filters the input data", {
  inputTab <- data.table(
    Sample1.PG.Quantity = c("1000", "Filtered", "2000", "0"),
    PG.ProteinGroups = c("P1", "P2", "P3", "P4"),
    PG.Genes = c("Gene1", "Gene2", "Gene3", "Gene4")
  )
  
  result <- readOneProteomDIA(inputTab, "Sample1")
  
  expect_equal(nrow(result), 2)
  expect_true(all(result$Intensity > 0))
})

test_that("readOneProteomDIA handles 'Filtered' values correctly", {
  inputTab <- data.table(
    Sample1.PG.Quantity = c("Filtered", "Filtered", "2000"),
    PG.ProteinGroups = c("P1", "P2", "P3"),
    PG.Genes = c("Gene1", "Gene2", "Gene3")
  )
  
  result <- readOneProteomDIA(inputTab, "Sample1")
  
  expect_equal(nrow(result), 1)
  expect_equal(result$Intensity, 2000)
})

test_that("readOneProteomDIA converts character values to numeric", {
  inputTab <- data.table(
    Sample1.PG.Quantity = c("1000", "Filtered", "2,456"),
    PG.ProteinGroups = c("P1", "P2", "P3"),
    PG.Genes = c("Gene1", "Gene2", "Gene3")
  )
  
  result <- readOneProteomDIA(inputTab, "Sample1")
  
  expect_equal(result$Intensity, c(1000, 2.456))
})

test_that("readOneProteomDIA renames columns correctly", {
  inputTab <- data.table(
    Sample1.PG.Quantity = c("1000", "Filtered", "2000"),
    PG.ProteinGroups = c("P1", "P2", "P3"),
    PG.Genes = c("Gene1", "Gene2", "Gene3")
  )
  
  result <- readOneProteomDIA(inputTab, "Sample1")
  
  expect_true(all(colnames(result) == c("Intensity", "UniprotID", "Gene")))
})

test_that("readOneProteomDIA returns NULL and issues a warning when all rows are filtered out", {
  inputTab <- data.table(
    Sample1.PG.Quantity = c("Filtered", "Filtered", "0"),
    PG.ProteinGroups = c("P1", "P2", "P3"),
    PG.Genes = c("Gene1", "Gene2", "Gene3")
  )
  
  expect_warning(result <- readOneProteomDIA(inputTab, "Sample1"))
  expect_null(result)
})

test_that("readOneProteomDIA ensures unique identifiers with no duplicates", {
  inputTab <- data.table(
    Sample1.PG.Quantity = c("1000", "Filtered", "2000"),
    PG.ProteinGroups = c("P1", "P2", "P3"),
    PG.Genes = c("Gene1", "Gene2", "Gene3")
  )
  
  result <- readOneProteomDIA(inputTab, "Sample1")
  
  expect_true(all(!duplicated(result$UniprotID)))
})



file <- system.file("extdata", "proteomeDIA_1.xls", package = "SmartPhos")

fileTable <- data.frame(type = "proteome", fileName = file, id = c("sample1", "sample2"))
# Test cases
test_that("readProteomeExperimentDIA returns a SummarizedExperiment", {
  
  result <- readProteomeExperimentDIA(fileTable)
  expect_s4_class(result, "SummarizedExperiment")
})

test_that("readProteomeExperimentDIA correctly filters and processes the input data", {
  
  result <- readProteomeExperimentDIA(fileTable)
  
  expect_equal(assay(result, "Intensity")["p1", "sample1"], 1000)
  expect_equal(assay(result, "Intensity")["p2", "sample2"], 2500)
  expect_equal(colnames(result), c("sample1", "sample2"))
  
})


test_that("readProteomeExperimentDIA returns NULL when no proteome data is present", {
  fileTable <- data.frame(type = "phosphoproteome", mfileName = file, id = "sample1")

  result <- readProteomeExperimentDIA(fileTable)

  expect_null(result)
})

test_that("readProteomeExperimentDIA throws an error when no proteins pass the threshold", {
  
  file <- system.file("extdata", "proteomeDIA_2.xls", package = "SmartPhos")
  
  fileTable <- data.frame(type = "proteome", fileName = file, id = "sample1")

  expect_error(readProteomeExperimentDIA(fileTable),
               "No proteins could pass the specified threshold in any sample!")
})

