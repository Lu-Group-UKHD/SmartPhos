########################### Tests for runFisher() ##############################

test_that("runFisher returns correct structure for non-PTM gene sets", {
  # Define genes of interest
  genes <- c("GeneA", "GeneB", "GeneC", "GeneD")

  # Define reference genes
  reference <- c("GeneE", "GeneF", "GeneG", "GeneH", "GeneI",
                 "GeneJ", "GeneK", "GeneL")

  # Define inputSet with gene set collections
  inputSet <- list(
    gsc = list(
      Set1 = c("GeneA", "GeneE", "GeneI"),
      Set2 = c("GeneB", "GeneF"),
      Set3 = c("GeneC", "GeneG", "GeneH", "GeneL")
    )
  )

  # Run the function
  result <- runFisher(genes, reference, inputSet, ptm = FALSE)

  # Check that result is a data frame
  expect_s3_class(result, "data.frame")

  # Check that result has the expected columns
  expect_true(all(c("Name", "Gene.number", "Set.size",
                    "pval", "padj", "Genes") %in% colnames(result)))

  # Check that 'Name' column matches the gene set names
  expect_equal(sort(result$Name), sort(names(inputSet$gsc)))

  # Expected 'Gene.number' per set
  expected_Gene_number <- c(1, 1, 1)

  expect_equal(result$Gene.number[order(result$Name)], expected_Gene_number)

  # Expected 'Set.size' per set
  expected_Set_size <- c(3, 2, 4)

  expect_equal(result$Set.size[order(result$Name)], expected_Set_size)

  # Expected intersecting genes per set
  expected_Genes <- list(Set1 = "GeneA", Set2 = "GeneB", Set3 = "GeneC")

  genSet <- result$Genes[order(result$Name)]
  names(genSet) <- names(inputSet$gsc)

  expect_equal(genSet, expected_Genes)

  # Check that 'pval' and 'padj' are numeric and within [0,1]
  expect_true(all(result$pval >= 0 & result$pval <= 1))
  expect_true(all(result$padj >= 0 & result$padj <= 1))
})

test_that("runFisher returns correct structure for PTM gene sets", {
  # Define genes of interest
  genes <- c("GeneA_pS1", "GeneB_pT2", "GeneC_pY3", "GeneD_pS4", "GeneK_pS11")

  # Define reference genes
  reference <- c("GeneE_pS5", "GeneF_pT6", "GeneG_pY7", "GeneH_pS8",
                 "GeneI_pT9", "GeneJ_pY10", "GeneK_pS11", "GeneL_pT12")

  # Define inputSet as a data frame with necessary columns
  inputSet <- data.frame(
    category = c(rep("CATEGORY1", 25), rep("KINASE", 5)), # 30 rows
    site.ptm = c(rep("p", 25), rep("other", 5)),
    signature = c(rep("Sig1", 10), rep("Sig2", 10), rep("Sig3", 10)),
    site.direction = c(rep("u", 5), rep("d", 5), rep("u", 5), rep("d", 5),
                       rep("u", 5), rep("d", 5)),
    site.annotation = paste0(
      c("GeneA_pS1", "GeneB_pS2", "GeneC_pY3", "GeneD_pS5", "GeneE_pS5",
        "GeneF_pT7", "GeneG_pY8", "GeneH_pS9", "GeneI_pT10", "GeneJ_pY11",
        "GeneK_pS1", "GeneL_pT12", "GeneM_pS13", "GeneD_pS4", "GeneO_pY3",
        "GeneP_pS1", "GeneQ_pT17", "GeneR_pY18", "GeneS_pS19", "GeneT_pT20",
        "GeneA_pS1", "GeneB_pS2", "GeneK_pS11", "GeneX_pY3", "GeneD_pS4",
        "GeneA_pS1", "GeneAA_pY27", "GeneAB_pS28", "GeneAC_pT29", "GeneAD_pY30"
      ),
      ":PMID123456"
    )
  )

  # Run the function
  result <- runFisher(genes, reference, inputSet, ptm = TRUE)

  # Check that result is a data frame
  expect_s3_class(result, "data.frame")

  # Check that result has the expected columns
  expect_true(all(c("Name", "Gene.number", "Set.size",
                    "pval", "padj", "Genes") %in% colnames(result)))

  # Expected signature names after processing
  expected_Names <- c("Sig1_upregulated", "Sig2_upregulated",
                      "Sig3_upregulated" )

  # Check that 'Name' column matches the expected Names
  expect_equal(sort(result$Name), sort(expected_Names))

  # Check that 'Gene.number' and 'Set.size' are non-negative integers
  expect_true(all(result$Gene.number >= 0))
  expect_true(all(result$Set.size >= result$Gene.number))

  # Check that 'pval' and 'padj' are numeric and within [0,1]
  expect_true(all(result$pval >= 0 & result$pval <= 1))
  expect_true(all(result$padj >= 0 & result$padj <= 1))

  # Check that 'Genes' column is a list
  expect_true(is.list(result$Genes))
})

test_that("runFisher throws an error when all gene sets are 'KINASE'", {
  # Define genes of interest
  genes <- c("GeneA_pS1", "GeneB_pT2", "GeneC_pY3", "GeneD_pS4")

  # Define reference genes
  reference <- c("GeneE_pS5", "GeneF_pT6", "GeneG_pY7", "GeneH_pS8",
                 "GeneI_pT9", "GeneJ_pY10", "GeneK_pS11", "GeneL_pT12")

  # Define inputSet as a data frame with 'KINASE' category only
  inputSet <- data.frame(
    category = rep("KINASE", 10),
    site.ptm = rep("p", 10),
    signature = rep("Sig1", 10),
    site.direction = c(rep("u", 5), rep("d", 5)),
    site.annotation = paste0(
      c("GeneA_pS1", "GeneB_pT2", "GeneC_pY3", "GeneD_pS4", "GeneE_pS5",
        "GeneF_pT6", "GeneG_pY7", "GeneH_pS8", "GeneI_pT9", "GeneJ_pY10"),
      ":PMID123456"
    )
  )

  # Since all 'category' is 'KINASE', function should throw an error
  expect_error(runFisher(genes, reference, inputSet, ptm = TRUE),
               "Genesets empty. Check inputSet argument.")
})

test_that("runFisher handles empty 'genes' vector", {
  # Define empty genes of interest
  genes <- character(0)

  # Define reference genes
  reference <- c("GeneE", "GeneF", "GeneG", "GeneH", "GeneI", "GeneJ",
                 "GeneK", "GeneL")

  # Define inputSet with gene set collections
  inputSet <- list(
    gsc = list(
      Set1 = c("GeneA", "GeneE", "GeneI"),
      Set2 = c("GeneB", "GeneF"),
      Set3 = c("GeneC", "GeneG", "GeneH", "GeneL")
    )
  )

  # Run the function
  result <- runFisher(genes, reference, inputSet, ptm = FALSE)

  # Since 'genes' is empty, 'Gene.number' should be zero for all sets
  expect_true(all(result$Gene.number == 0))

  # p-values should be 1, as no enrichment
  expect_true(all(result$pval == 1))

  # Adjusted p-values should also be 1
  expect_true(all(result$padj == 1))

  # 'Genes' column should be empty lists
  expect_true(all(sapply(result$Genes, length) == 0))
})

test_that("runFisher handles 'genes' equal to 'reference'", {
  # Define genes of interest
  genes <- c("GeneA", "GeneB", "GeneC", "GeneD")

  # Define reference genes as same as 'genes'
  reference <- c("GeneA", "GeneB", "GeneC", "GeneD")

  # Define inputSet with gene set collections
  inputSet <- list(
    gsc = list(
      Set1 = c("GeneA", "GeneD", "GeneI"),
      Set2 = c("GeneB", "GeneF"),
      Set3 = c("GeneC", "GeneG", "GeneH", "GeneL")
    )
  )

  # Run the function
  result <- runFisher(genes, reference, inputSet, ptm = FALSE)

  # p-values should be 1, as no enrichment
  expect_true(all(result$pval == 1))

  # Adjusted p-values should also be 1
  expect_true(all(result$padj == 1))

  # 'Gene.number' should be correct
  expected_Gene_number <- c(2, 1, 1)

  expect_equal(result$Gene.number[order(result$Name)], expected_Gene_number)
})

test_that("runFisher correctly adjusts p-values", {
  # Define genes of interest
  genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")

  # Define reference genes
  reference <- c("GeneF", "GeneG", "GeneH", "GeneI", "GeneJ", "GeneK",
                 "GeneL", "GeneM", "GeneN", "GeneO", "GeneP")

  # Define inputSet with gene set collections
  inputSet <- list(
    gsc = list(
      Set1 = c("GeneA", "GeneF", "GeneI"),
      Set2 = c("GeneB", "GeneG", "GeneJ"),
      Set3 = c("GeneC", "GeneH", "GeneK"),
      Set4 = c("GeneD", "GeneL", "GeneM"),
      Set5 = c("GeneE", "GeneN", "GeneO"),
      Set6 = c("GeneX", "GeneY", "GeneZ") # No overlap
    )
  )

  # Run the function
  result <- runFisher(genes, reference, inputSet, ptm = FALSE)

  # Check that 'padj' is equal to p.adjust(pval, method = "BH")
  expected_padj <- p.adjust(result$pval, method = "BH")

  expect_equal(result$padj, expected_padj)
})

test_that("runFisher correctly identifies intersecting genes", {
  # Define genes of interest
  genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")

  # Define reference genes
  reference <- c("GeneF", "GeneG", "GeneH", "GeneI", "GeneJ", "GeneK",
                 "GeneL", "GeneM", "GeneN", "GeneO", "GeneP")

  # Define inputSet with gene set collections
  inputSet <- list(
    gsc = list(
      Set1 = c("GeneA", "GeneF", "GeneI"),
      Set2 = c("GeneB", "GeneG", "GeneJ"),
      Set3 = c("GeneC", "GeneH", "GeneK"),
      Set4 = c("GeneD", "GeneL", "GeneM"),
      Set5 = c("GeneE", "GeneN", "GeneO"),
      Set6 = c("GeneX", "GeneY", "GeneZ") # No overlap
    )
  )

  # Run the function
  result <- runFisher(genes, reference, inputSet, ptm = FALSE)

  # Expected 'Genes' column
  expected_Genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")

  expect_equal(unlist(result$Genes[order(result$Name)]), expected_Genes)
})

test_that("runFisher filters out gene sets with Set.size == 0", {
  # Define genes of interest
  genes <- c("GeneA", "GeneB", "GeneC", "GeneD", "GeneE")

  # Define reference genes
  reference <- c("GeneF", "GeneG", "GeneH", "GeneI", "GeneJ", "GeneK",
                 "GeneL", "GeneM", "GeneN", "GeneO", "GeneP")

  # Define inputSet with gene set collections
  inputSet <- list(
    gsc = list(
      Set1 = c("GeneA", "GeneF", "GeneI"),
      Set2 = c("GeneB", "GeneG", "GeneJ"),
      Set3 = c("GeneC", "GeneH", "GeneK"),
      Set4 = c("GeneD", "GeneL", "GeneM"),
      Set5 = c("GeneE", "GeneN", "GeneO"),
      Set6 = c("GeneX", "GeneY", "GeneZ") # No overlap
    )
  )

  # Run the function
  result <- runFisher(genes, reference, inputSet, ptm = FALSE)

  # 'Set6' has no genes in 'genes' or 'reference', so 'Set.size' == 0, and
  # should be filtered out
  expect_false("Set6" %in% result$Name)

  # Other sets should be present
  expect_true(all(c("Set1", "Set2", "Set3", "Set4", "Set5") %in% result$Name))
})


####################### Tests for enrichDifferential() #########################

library(piano)
# Load multiAssayExperiment object
data("dia_example")
# Get SummarizedExperiment object
se <- dia_example[["Phosphoproteome"]]
colData(se) <- colData(dia_example)
newSE <- preprocessPhos(se, normalize = TRUE)
# Perform differtial expression analyis
dea <- performDifferentialExp(se = newSE, assay = "Intensity",
                              method = "limma", reference = "1stCrtl",
                              target = "EGF",condition = "treatment")
# Load gene set
genesetPath <- appDir <- system.file("shiny-app/geneset",package = "SmartPhos")
inGMT <- loadGSC(paste0(genesetPath,"/Cancer_Hallmark.gmt"),type="gmt")


test_that("enrichDifferential handles pathway enrichment with PAGE", {

    result <- enrichDifferential(dea = dea$resDE, type = "Pathway enrichment",
                                 gsaMethod = "PAGE", geneSet = inGMT,
                                 statType = "stat", nPerm = 200,
                                 sigLevel = 0.05, ifFDR = FALSE)

    # Test output structure
    expect_s3_class(result, "data.frame")
    expect_true(all(c("Name", "Stat", "p.up", "p.down") %in% colnames(result)))
})

test_that("enrichDifferential respects significance threshold", {

    result <- enrichDifferential(dea = dea$resDE, type = "Pathway enrichment",
                                 gsaMethod = "PAGE", geneSet = inGMT,
                                 statType = "stat", nPerm = 200,
                                 sigLevel = 0.01, ifFDR = FALSE)

    expect_true(all(result$p.up <= 0.01 | result$p.down <= 0.01))
})


######################### Tests for clusterEnrich() ############################

# Mock data for testing
set.seed(123)
mock_clusterTab <- data.frame(
  feature = paste0("gene", 1:10),
  cluster = c(rep("Cluster1", 5), rep("Cluster2", 5))
)

mock_se <- function() {
  # Create assay data
  assay_data <- matrix(abs(rnorm(1000, 10, 2)), nrow = 100, ncol = 10)
  colnames(assay_data) <- paste0("Sample", 1:10)
  rownames(assay_data) <- paste0("gene", 1:100)

  # Create sample metadata
  sample_data <- data.frame(sample = colnames(assay_data),
                            sampleType = rep(c("FullProteome", "Phospho"), each = 5),
                            group = rep(c("A", "B"), each = 5))
  rownames(sample_data) <- colnames(assay_data)

  # Create gene metadata
  gene_data <- data.frame(UniprotID = rownames(assay_data),
                          Gene = paste0("gene", 1:100),
                          site = paste0("site", 1:10))
  rownames(gene_data) <- rownames(assay_data)

  # Create SummarizedExperiment object
  se <- SummarizedExperiment(assays = list(Intensity = assay_data),
                             colData = sample_data, rowData = gene_data)
  return(se)
}

mock_inputSet <- list(
  gsc = list(
    pathway1 = c("gene1", "gene2", "gene3"),
    pathway2 = c("gene6", "gene7", "gene8")
  )
)

mock_reference <- c("gene1", "gene2", "gene3", "gene4", "gene5",
                    "gene6", "gene7", "gene8", "gene9", "gene10")

# Test that the function runs without errors with default parameters
test_that("clusterEnrich runs without errors with default parameters", {
  result <- clusterEnrich(
    clusterTab = mock_clusterTab,
    se = mock_se(),
    inputSet = mock_inputSet
  )

  expect_true(is.list(result))
  expect_true("table" %in% names(result))
  expect_true("plot" %in% names(result))
  expect_s3_class(result$plot, "ggplot")
})

# Test that the correct number of clusters is returned
test_that("clusterEnrich returns correct number of clusters", {
  result <- clusterEnrich(
    clusterTab = mock_clusterTab,
    se = mock_se(),
    inputSet = mock_inputSet
  )

  unique_clusters <- unique(mock_clusterTab$cluster)
  result_clusters <- unique(result$table$cluster)

  expect_equal(length(unique_clusters), length(result_clusters))
  expect_equal(sort(unique_clusters), sort(result_clusters))
})

# Test that the function correctly filters results based on p-value threshold
test_that("clusterEnrich filters results based on p-value threshold", {
  result <- clusterEnrich(
    clusterTab = mock_clusterTab,
    se = mock_se(),
    inputSet = mock_inputSet,
    filterP = 0.01
  )

  expect_true(all(result$table$pval <= 0.01))
})

# Test that FDR adjustment is correctly applied
test_that("clusterEnrich applies FDR adjustment correctly", {
  result <- clusterEnrich(
    clusterTab = mock_clusterTab,
    se = mock_se(),
    inputSet = mock_inputSet,
    ifFDR = TRUE,
    filterP = 0.05
  )

  expect_true(all(result$table$padj <= 0.05))
})


# Test PTM-specific gene sets analysis
test_that("clusterEnrich handles PTM-specific gene sets correctly", {

  inputSet <- data.frame(
    category = c(rep("CATEGORY1", 25), rep("KINASE", 5)), # 30 rows
    site.ptm = c(rep("p", 25), rep("other", 5)),
    signature = c(rep("Sig1", 5), rep("Sig2", 10), rep("Sig3", 15)),
    site.direction = c(rep("u", 5), rep("d", 5), rep("u", 5), rep("d", 5),
                       rep("u", 5), rep("d", 5)),
    site.annotation = paste0("site", 1:30)
    )

  result <- clusterEnrich(
    clusterTab = mock_clusterTab,
    se = mock_se(),
    inputSet = inputSet,
    ptm = TRUE
  )

  expect_true(is.list(result))
  expect_true("table" %in% names(result))
  expect_true("plot" %in% names(result))
  expect_s3_class(result$plot, "ggplot")
})
