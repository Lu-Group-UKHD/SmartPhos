###################### Tests for readOnePhosDIA() #############################

# Mock data for testing
inputTab <- data.table(
  Sample_1.PTM.SiteProbability = c("0.8", "Filtered", "0.9", "0.6"),
  Sample_1.PTM.Quantity = c("1000", "Filtered", "2000", "3000"),
  PTM.CollapseKey = c("Key1", "Key2", "Key3", "Key4"),
  PG.UniProtIds = c("P1", "P2", "P3", "P4"),
  PG.Genes = c("Gene1", "Gene2", "Gene3", "Gene4"),
  PTM.Multiplicity = c(1, 1, 2, 1),
  PTM.SiteLocation = c("Loc1", "Loc2", "Loc3", "Loc4"),
  PTM.SiteAA = c("A", "B", "C", "D"),
  PTM.FlankingRegion = c("Region1", "Region2", "Region3", "Region4")
)

# Unit tests
test_that("readOnePhosDIA processes data correctly", {
  result <- readOnePhosDIA(inputTab, "Sample_1", localProbCut = 0.75, 
                           removeDup = FALSE)
  expect_equal(nrow(result), 2)
  expect_true(all(result$Intensity > 0))
  expect_true(all(result$PG.UniProtIds %in% c("P1", "P3")))
})

test_that("readOnePhosDIA handles no rows passing filters", {
  expect_warning(result <- readOnePhosDIA(inputTab, "Sample_1", 
                                          localProbCut = 0.95, 
                                          removeDup = FALSE))
  expect_null(result)
})

test_that("readOnePhosDIA handles removeDup correctly", {
  inputTab_dup <- rbind(inputTab, inputTab[1,])
  inputTab_dup$PG.UniProtIds[5] <- "P1"
  inputTab_dup$PTM.Quantity.Sample_1[5] <- "1000"
  result <- readOnePhosDIA(inputTab_dup, "Sample_1", localProbCut = 0.75, 
                           removeDup = TRUE)
  expect_equal(nrow(result), 2)
  expect_true(all(result$PG.UniProtIds %in% c("P1", "P3")))
})

test_that("readOnePhosDIA checks for required columns", {
  inputTab_missing <- inputTab[, -"Sample_1.PTM.SiteProbability", with = FALSE]
  expect_error(readOnePhosDIA(inputTab_missing, "Sample_1"), 
               "Sample not found in quantification file")
})

test_that("readOnePhosDIA handles multiplicity correctly", {
  inputTab <- data.table(
    Sample_1.PTM.SiteProbability = c("0.8", "Filtered", "0.9", "0.85"),
    Sample_1.PTM.Quantity = c("1000", "Filtered", "2000", "1500"),
    PTM.CollapseKey = c("P1_S23_M1", "P2_S56_M1", "P1_S23_M2", "P4_S98_M1"),
    PG.UniProtIds = c("P1", "P2", "P1", "P4"),
    PG.Genes = c("Gene1", "Gene2", "Gene3", "Gene4"),
    PTM.Multiplicity = c(1, 1, 2, 1),
    PTM.SiteLocation = c("Loc1", "Loc2", "Loc3", "Loc4"),
    PTM.SiteAA = c("A", "B", "C", "D"),
    PTM.FlankingRegion = c("Region1", "Region2", "Region3", "Region4")
  )

  result <- readOnePhosDIA(inputTab, "Sample_1", localProbCut = 0.75, 
                           removeDup = FALSE)

  expect_equal(result$Intensity[1], 3000)
  expect_equal(result$Intensity[2], 1500)
})

test_that("readOnePhosDIA handles multiple samples correctly", {
  inputTab_multi <- data.table(
    Sample_1.PTM.SiteProbability = c("0.8", "Filtered", "0.9", "0.6"),
    Sample_1.PTM.Quantity = c("1000", "Filtered", "2000", "3000"),
    Sample_2.PTM.SiteProbability = c("0.7", "0.8", "Filtered", "0.9"),
    Sample_2.PTM.Quantity = c("1500", "2500", "Filtered", "3500"),
    PTM.CollapseKey = c("Key1", "Key2", "Key3", "Key4"),
    PG.UniProtIds = c("P1", "P2", "P3", "P4"),
    PG.Genes = c("Gene1", "Gene2", "Gene3", "Gene4"),
    PTM.Multiplicity = c(1, 1, 2, 1),
    PTM.SiteLocation = c("Loc1", "Loc2", "Loc3", "Loc4"),
    PTM.SiteAA = c("A", "B", "C", "D"),
    PTM.FlankingRegion = c("Region1", "Region2", "Region3", "Region4")
  )

  result1 <- readOnePhosDIA(inputTab_multi, "Sample_1", localProbCut = 0.75, 
                            removeDup = FALSE)
  result2 <- readOnePhosDIA(inputTab_multi, "Sample_2", localProbCut = 0.75, 
                            removeDup = FALSE)

  expect_equal(nrow(result1), 2)
  expect_equal(nrow(result2), 2)
  expect_true(all(result1$Intensity > 0))
  expect_true(all(result2$Intensity > 0))
})


#################### Tests for readPhosphoExperimentDIA() ######################

file <- system.file("extdata", "phosDIA_1.xls", package = "SmartPhos")
fileTable <- data.frame(type = "phosphoproteome", fileName = file, 
                        id = c("Sample_1"))

test_that("readPhosphoExperimentDIA processes data correctly", {
  
  result <- readPhosphoExperimentDIA(fileTable, localProbCut = 0.75, 
                                     onlyReviewed = FALSE, 
                                     showProgressBar = FALSE)
  
  expect_s4_class(result, "SummarizedExperiment")
  expect_equal(ncol(assay(result)), 1)
  expect_equal(nrow(assay(result)), 2)
  expect_true(all(rowData(result)$UniprotID %in% c("P1", "P3")))
})

test_that("readPhosphoExperimentDIA handles no rows passing filters", {

  expect_error(readPhosphoExperimentDIA(fileTable, localProbCut = 0.95, 
                                        onlyReviewed = FALSE, 
                                        showProgressBar = FALSE),
    "No phosphorylation site could pass the specified threshold in any sample!"
  )
})

file <- system.file("extdata", "phosDIA_2.xls", package = "SmartPhos")
fileTable <- data.frame(type = "phosphoproteome", fileName = file, id = c("Sample_1"))

test_that("readPhosphoExperimentDIA handles missing PG.ProteinGroups and PG.UniProtIds", {
  
  expect_error(readPhosphoExperimentDIA(fileTable, localProbCut = 0.75, 
                                        onlyReviewed = FALSE, 
                                        showProgressBar = FALSE),
               "Either PG.ProteinGroups or PG.UniProtIds should be in the quantification table"
  )
})

###################### Tests for readOneProteomeDIA() ##########################

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

test_that("readOneProteomDIA checks for required columns", {
  inputTab <- data.table(
    Sample1.PG.Quantity = c("1000", "Filtered", "2000"),
    PG.ProteinGroups = c("P1", "P2", "P3"),
    PG.Genes = c("Gene1", "Gene2", "Gene3")
  )
  expect_error(readOnePhosDIA(inputTab, "Sample_2"), "Sample not found in quantification file")
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


################### Tests for readProteomeExperimentDIA() ######################


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

