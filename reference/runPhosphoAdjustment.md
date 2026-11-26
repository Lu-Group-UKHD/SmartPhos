# Run Phospho Adjustment

`runPhosphoAdjustment` performs phospho adjustment on a
`MultiAssayExperiment` object to normalize the phosphoproteome data.

## Usage

``` r
runPhosphoAdjustment(
  maeData,
  normalization = FALSE,
  minOverlap = 3,
  completeness = 0,
  ncore = 1
)
```

## Arguments

- maeData:

  A `MultiAssayExperiment` object containing phosphoproteome and
  proteome data.

- normalization:

  A `logical` value indicating whether to perform normalization. Default
  is `FALSE`.

- minOverlap:

  A `numeric` value specifying the minimum number of overlapping
  peptides required between samples. Default is 3.

- completeness:

  A `numeric` value indicating the required completeness of data for
  features to be included. Default is 0.

- ncore:

  A `numeric` value specifying the number of cores to use for parallel
  processing. Default is 1.

## Value

A `MultiAssayExperiment` object with adjusted phosphoproteome data.

## Details

The function performs the following steps:

1.  Defines an optimization function to minimize the sum of squared
    differences between pairs of samples.

2.  Calculates the ratio matrix of phosphoproteome to total/full
    proteome data.

3.  Subsets features based on completeness criteria.

4.  Performs a sanity check to identify and exclude problematic samples.

5.  Sets initial values for the adjustment factor based on column
    medians.

6.  Estimates the adjustment factor using parallel optimization.

7.  Adjusts the phosphoproteome measurements using the estimated
    adjustment factor.
