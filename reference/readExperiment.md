# Read and Process the DDA experiment.

`readExperiment` reads and processes DDA (Data-Dependent Acquisition)
phosphoproteomic and proteomic data from a given file table, and returns
a `MultiAssayExperiment` object.

## Usage

``` r
readExperiment(
  fileTable,
  localProbCut = 0.75,
  scoreDiffCut = 5,
  fdrCut = 0.1,
  scoreCut = 10,
  pepNumCut = 1,
  ifLFQ = TRUE,
  annotation_col = c(),
  verbose = FALSE
)
```

## Arguments

- fileTable:

  A `data.frame` containing information about the input files, including
  searchType, id, sample, and other annotations.

- localProbCut:

  `Numeric`, local probability cutoff for filtering phosphoproteomic
  data. Default is 0.75.

- scoreDiffCut:

  `Numeric`, score difference cutoff for filtering phosphoproteomic
  data. Default is 5.

- fdrCut:

  `Numeric`, false discovery rate cutoff for filtering proteomic data.
  Default is 0.1.

- scoreCut:

  `Numeric`, score cutoff for filtering proteomic data. Default is 10.

- pepNumCut:

  `Numeric`, peptide number cutoff for filtering proteomic data. Default
  is 1.

- ifLFQ:

  `Logical`, whether to use LFQ quantification for proteomic data.
  Default is `TRUE`.

- annotation_col:

  A `character` vector specifying additional columns to be included in
  the sample annotation. Default is an empty vector.

- verbose:

  `Logical` value indicating whether to print detailed information.
  Default is `FALSE`.

## Value

A `MultiAssayExperiment` object containing the processed
phosphoproteomic and proteomic data from a DDA experiment.

## Details

The function performs the following steps:

- Reads and processes the phosphoproteomic data using the
  `readPhosphoExperiment` function.

- Reads and processes the proteomic data using the
  `readProteomeExperiment` function.

- Prepares the sample annotation table.

- Constructs and returns a `MultiAssayExperiment` object containing the
  processed data.

## Examples

``` r
# Example usage:
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
readExperiment(fileTable, localProbCut = 0.75, scoreDiffCut = 5,
fdrCut = 0.1, scoreCut = 10, pepNumCut = 1, ifLFQ = TRUE,
annotation_col = c("id"))
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] Phosphoproteome: SummarizedExperiment with 2 rows and 1 columns
#>  [2] Proteome: SummarizedExperiment with 2 rows and 1 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
