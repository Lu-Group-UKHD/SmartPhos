# Read Phosphorylation Experiment Data from DIA

`readPhosphoExperimentDIA` reads and processes phosphorylation data from
DIA experiments, applying filters for localization probability, and
optionally including only reviewed proteins. It constructs a
`SummarizedExperiment` object.

## Usage

``` r
readPhosphoExperimentDIA(
  fileTable,
  localProbCut = 0.75,
  onlyReviewed = TRUE,
  showProgressBar = FALSE
)
```

## Arguments

- fileTable:

  A `data.table` or `data.frame` containing metadata about the files to
  read. Must include columns fileName, searchType, and optionally
  outputID.

- localProbCut:

  A `numeric` value specifying the cutoff for localization probability.
  Default is 0.75.

- onlyReviewed:

  A `logical` value indicating whether to include only reviewed
  proteins. Default is `TRUE`.

- showProgressBar:

  A `logical` value indicating whether to show a progress bar. Default
  is `FALSE`.

## Value

A `SummarizedExperiment` object containing the processed phosphorylation
data.

## Details

This function processes phosphorylation data from DIA experiments by
filtering based on localization probability and non-zero intensity,
handling multiplicity, and optionally including only reviewed proteins.
The resulting data is returned as a `SummarizedExperiment` object with
annotations and an intensity matrix.

## Examples

``` r
file <- system.file("extdata", "phosDIA_1.xls", package = "SmartPhos")
fileTable <- data.frame(searchType = "phosphoproteome", fileName = file,
id = c("Sample_1"))
readPhosphoExperimentDIA(fileTable, localProbCut = 0.75,
onlyReviewed = FALSE, showProgressBar = FALSE)
#> class: SummarizedExperiment 
#> dim: 2 1 
#> metadata(0):
#> assays(1): Intensity
#> rownames(2): s1 s2
#> rowData names(7): UniprotID Gene ... Sequence site
#> colnames(1): Sample_1
#> colData names(0):
```
