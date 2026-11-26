# Read and Process Proteomics Experiment Data

`readProteomeExperiment` reads and processes proteomics data from
multiple samples, applying various quality filters, and returns a
`SummarizedExperiment` object.

## Usage

``` r
readProteomeExperiment(
  fileTable,
  fdrCut = 0.1,
  scoreCut = 10,
  pepNumCut = 1,
  ifLFQ = TRUE
)
```

## Arguments

- fileTable:

  A `data.table` or `data.frame` containing the file information with
  columns for file names, sample names, and IDs.

- fdrCut:

  A `numeric` value specifying the maximum false discovery rate (FDR)
  threshold. Default is 0.1.

- scoreCut:

  A `numeric` value specifying the minimum score threshold. Default is
  10.

- pepNumCut:

  A `numeric` value specifying the minimum number of peptides required
  for a protein to be included. Default is 1.

- ifLFQ:

  A `logical` value indicating whether to use LFQ quantification.
  Default is `TRUE`.

## Value

A `SummarizedExperiment` object containing the processed proteomics
data.

## Details

This function processes proteomics data by filtering based on FDR,
score, and peptide count, and optionally using LFQ quantification. It
aggregates the data from multiple samples and constructs a
`SummarizedExperiment` object.

## Examples

``` r
file1 <- system.file("extdata", "phosDDA_1.xls", package = "SmartPhos")
file2 <- system.file("extdata", "proteomeDDA_1.xls", package = "SmartPhos")
# Create fileTable
fileTable <- data.frame(
   searchType = c("phosphoproteome", "proteome"),
   fileName = c(file1, file2),
   sample = c("Sample1", "sample1"),
   id = c("s1", "s2")
)
# Call the function
readProteomeExperiment(fileTable, fdrCut = 0.1, scoreCut = 10,
pepNumCut = 1, ifLFQ = TRUE)
#> class: SummarizedExperiment 
#> dim: 2 1 
#> metadata(0):
#> assays(1): Intensity
#> rownames(2): p1 p2
#> rowData names(3): UniprotID Gene PeptideCounts
#> colnames(1): s2
#> colData names(0):
```
