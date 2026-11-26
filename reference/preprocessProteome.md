# Preprocess Proteome Data

`preprocessProteome` preprocesses proteome data stored in a
`SummarizedExperiment` object by performing filtering, transformation,
normalization, imputation, and batch effect removal.

## Usage

``` r
preprocessProteome(
  seData,
  filterList = NULL,
  missCut = 50,
  transform = c("log2", "vst", "none"),
  normalize = FALSE,
  getPP = FALSE,
  removeOutlier = NULL,
  impute = c("none", "QRILC", "MLE", "bpca", "missForest", "MinDet"),
  batch = NULL,
  verbose = FALSE,
  scaleFactorTab = NULL
)
```

## Arguments

- seData:

  A `SummarizedExperiment` object containing proteome data.

- filterList:

  A `list` of filters to apply on the samples. Default is `NULL`.

- missCut:

  `Numeric` value specifying the missing value cutoff percentage for
  filtering features. Default is 50.

- transform:

  `Character` string specifying the transformation method ("log2",
  "vst", "none"). Default is "log2".

- normalize:

  `Logical` value indicating whether to normalize the data. Default is
  `FALSE`.

- getPP:

  `Logical` value indicating whether to retrieve PP samples. Default is
  `FALSE`.

- removeOutlier:

  `Character` vector of samples to be removed as outliers. Default is
  `NULL`.

- impute:

  `Character` string specifying the imputation method ("QRILC", "MLE",
  "bpca", "missForest", "MinDet", "none"). Default is "none".

- batch:

  `Character` vector specifying batch effects to remove. Default is
  `NULL`.

- verbose:

  `Logical` value indicating whether to print detailed information.
  Default is `FALSE`.

- scaleFactorTab:

  `Data frame` containing scale factors for normalization. Default is
  `NULL`.

## Value

A `SummarizedExperiment` object with preprocessed proteome data.

## Examples

``` r
library(SummarizedExperiment)
# Load multiAssayExperiment object
data("dia_example")
# Get SummarizedExperiment object
se <- dia_example[["Proteome"]]
colData(se) <- colData(dia_example)
# Call the function
preprocessProteome(seData = se, normalize = TRUE, impute = "QRILC")
#> Imputing along margin 2 (samples/columns).
#> class: SummarizedExperiment 
#> dim: 491 27 
#> metadata(0):
#> assays(2): Intensity imputed
#> rownames(491): p1 p2 ... p499 p500
#> rowData names(2): UniprotID Gene
#> colnames(27): FullProteome_1stCrtl_0min_rep2
#>   FullProteome_1stCrtl_0min_rep3 ... FullProteome_HGF_24h_rep1
#>   FullProteome_HGF_100min_rep1
#> colData names(6): sample treatment ... sampleType sampleName
```
