# Normalize Phosphoproteome by Total Proteome / Full Proteome

`normByFullProteome` normalizes the phosphoproteome data by the
corresponding total proteome / full proteome data in a
`MultiAssayExperiment` object. The "Phosphoproteome" assay in the
MultiAssayExperiment will be replaced by the ratio.

## Usage

``` r
normByFullProteome(mae, replace = TRUE)
```

## Arguments

- mae:

  A `MultiAssayExperiment` object containing both phosphoproteome and
  proteome assays.

- replace:

  `Logical`, whether to replace the existing phosphoproteome assay in
  the `MultiAssayExperiment` object. If replace = `FALSE`, a new assay,
  phosphoRatio, will be created. Default is `TRUE`.

## Value

A `MultiAssayExperiment` object with the normalized phosphoproteome
data.

## Details

The function performs the following steps:

- Checks if both phosphoproteome and proteome assays are present in the
  `MultiAssayExperiment` object.

- Extracts the phosphoproteome and proteome assays along with the sample
  annotations.

- Matches the samples between the phosphoproteome and proteome assays.

- Normalizes the phosphoproteome data by dividing it by the
  corresponding proteome data.

- Replaces the phosphoproteome assay in the `MultiAssayExperiment`
  object or adds the normalized data as a new assay, depending on the
  `replace` parameter.

## Examples

``` r
# load mae object
data("dia_example")
# call the function
normByFullProteome(dia_example)
#> harmonizing input:
#>   removing 27 sampleMap rows with 'colname' not in colnames of experiments
#> A MultiAssayExperiment object of 2 listed
#>  experiments with user-defined names and respective classes.
#>  Containing an ExperimentList class object of length 2:
#>  [1] Phosphoproteome: SummarizedExperiment with 4 rows and 20 columns
#>  [2] Proteome: SummarizedExperiment with 500 rows and 47 columns
#> Functionality:
#>  experiments() - obtain the ExperimentList instance
#>  colData() - the primary/phenotype DataFrame
#>  sampleMap() - the sample coordination DataFrame
#>  `$`, `[`, `[[` - extract colData columns, subset, or experiment
#>  *Format() - convert into a long or wide DataFrame
#>  assays() - convert ExperimentList to a SimpleList of matrices
#>  exportClass() - save data to flat files
```
