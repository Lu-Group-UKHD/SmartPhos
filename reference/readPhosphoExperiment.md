# Read Phosphorylation Experiment Data

`readPhosphoExperiment` reads and processes phosphorylation experiment
data from multiple files, filtering based on localization probability
and score difference, and constructs a `SummarizedExperiment` object.

## Usage

``` r
readPhosphoExperiment(fileTable, localProbCut = 0.75, scoreDiffCut = 5)
```

## Arguments

- fileTable:

  A `data.table` or `data.frame` containing information about the files,
  including columns for file names, sample names, and other relevant
  metadata. It must include a column named "searchType" with value
  "phosphoproteome" for relevant entries.

- localProbCut:

  A `numeric` value specifying the cutoff for localization probability.
  Default is 0.75.

- scoreDiffCut:

  A `numeric` value specifying the cutoff for score difference. Default
  is 5.

## Value

A `SummarizedExperiment` object containing the processed phosphorylation
data.

## Details

This function reads phosphorylation data from multiple files as
specified in fileTable, filters the data based on localization
probability and score difference, and removes reverse and potential
contaminant entries. It constructs an intensity matrix and annotation
data, which are then used to create a `SummarizedExperiment` object.

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
readPhosphoExperiment(fileTable, localProbCut = 0.75, scoreDiffCut = 5)
#> class: SummarizedExperiment 
#> dim: 2 1 
#> metadata(0):
#> assays(1): Intensity
#> rownames(2): s1 s2
#> rowData names(6): UniprotID Gene ... Sequence site
#> colnames(1): s1
#> colData names(0):
```
