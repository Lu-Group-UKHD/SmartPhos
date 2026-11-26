# Read and Process a DIA Experiment

`readExperimentDIA` reads and processes DIA (Data-Independent
Acquisition) data for both phosphoproteome and proteome experiments, and
constructs a `MultiAssayExperiment` object.

## Usage

``` r
readExperimentDIA(
  fileTable,
  localProbCut = 0.75,
  annotation_col = c(),
  onlyReviewed = TRUE,
  normalizeByProtein = FALSE,
  verbose = FALSE
)
```

## Arguments

- fileTable:

  A `data frame` containing metadata about the files to be read. Must
  contain columns searchType, fileName, id, and optionally outputID.

- localProbCut:

  `Numeric`, the local probability cutoff for phosphoproteomic data.
  Default is 0.75.

- annotation_col:

  A `character` vector specifying the columns in fileTable to be used
  for sample annotation.

- onlyReviewed:

  A `logical` value indicating whether to include only reviewed
  proteins. Default is `TRUE`.

- normalizeByProtein:

  `Logical`, whether to normalize the data by protein. Default is
  `FALSE`.

- verbose:

  `Logical` value indicating whether to print detailed information.
  Default is `FALSE`.

## Value

A `MultiAssayExperiment` object containing the processed phosphoproteome
and proteome data.

## Details

The function performs the following steps:

- Reads and processes phosphoproteomic data using
  `readPhosphoExperimentDIA`.

- Reads and processes proteomic data using `readProteomeExperimentDIA`.

- Prepares sample annotations based on the provided fileTable and
  annotation_col.

- Constructs a `MultiAssayExperiment` object with the processed data and
  sample annotations.

The `readPhosphoExperimentDIA` and `readProteomeExperimentDIA` functions
are used to read and filter the data for phosphoproteome and proteome
experiments, respectively, and they must be available in the
environment.

## Examples

``` r
# Example usage:
file1 <- system.file("extdata", "phosDIA_1.xls", package = "SmartPhos")
file2 <- system.file("extdata", "proteomeDIA_1.xls", package = "SmartPhos")
# Create fileTable
fileTable <- data.frame(
  searchType = c("phosphoproteome", "proteome", "proteome"),
  fileName = c(file1, file2, file2),
  id = c("Sample_1", "sample1", "sample2"),
  outputID = c("s1", "s2", "s3")
)
# Call the function
readExperimentDIA(fileTable, localProbCut = 0.75, annotation_col = c("id"),
onlyReviewed = FALSE, normalizeByProtein = FALSE)
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] Phosphoproteome: SummarizedExperiment with 2 rows and 1 columns
#>  [2] Proteome: SummarizedExperiment with 2 rows and 2 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
