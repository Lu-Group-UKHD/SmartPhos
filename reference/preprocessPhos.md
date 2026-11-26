# Preprocess Phosphoproteome Data

`preprocessPhos` preprocesses phosphoproteome data stored in a
`SummarizedExperiment` object by performing filtering, transformation,
normalization, imputation, and batch effect removal.

## Usage

``` r
preprocessPhos(
  seData,
  filterList = NULL,
  missCut = 50,
  transform = c("log2", "vst", "none"),
  normalize = FALSE,
  getFP = FALSE,
  removeOutlier = NULL,
  assayName = NULL,
  batch = NULL,
  scaleFactorTab = NULL,
  impute = c("none", "QRILC", "MLE", "bpca", "missForest", "MinDet"),
  verbose = FALSE
)
```

## Arguments

- seData:

  A `SummarizedExperiment` object containing phosphoproteome data.

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

- getFP:

  `Logical` value indicating whether to retrieve FP samples. Default is
  `FALSE`.

- removeOutlier:

  `Character` vector of samples to be removed as outliers. Default is
  `NULL`.

- assayName:

  `Character` string specifying the assay name in the
  SummarizedExperiment object. Default is `NULL`.

- batch:

  `Character` vector specifying batch effects to remove. Default is
  `NULL`.

- scaleFactorTab:

  `Data frame` containing scale factors for normalization. Default is
  `NULL`.

- impute:

  `Character` string specifying the imputation method ("QRILC", "MLE",
  "bpca", "missForest", "MinDet", "none"). Default is "none".

- verbose:

  `Logical` value indicating whether to print detailed information.
  Default is `FALSE`.

## Value

A `SummarizedExperiment` object with preprocessed phosphoproteome data.

## Examples

``` r
library(SummarizedExperiment)
# Load multiAssayExperiment object
data("dia_example")
# Get SummarizedExperiment object
se <- dia_example[["Phosphoproteome"]]
colData(se) <- colData(dia_example)
# Call the function
preprocessPhos(seData = se, normalize = TRUE, impute = "QRILC")
#> Imputing along margin 2 (samples/columns).
#> class: SummarizedExperiment 
#> dim: 172 20 
#> metadata(0):
#> assays(2): Intensity imputed
#> rownames(172): s4 s6 ... s499 s500
#> rowData names(7): UniprotID Gene ... Sequence site
#> colnames(20): Phospho_1stCrtl_0min_rep2 Phospho_1stCrtl_0min_rep3 ...
#>   Phospho_HGF_24h_rep1 Phospho_HGF_100min_rep1
#> colData names(6): sample treatment ... sampleType sampleName
```
