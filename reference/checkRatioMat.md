# Check the PP/FP ratio matrix and remove feature that do not meet requirements

`checkRatioMat` checks the ratio `matrix` for samples that do not have
sufficient overlap of phospho-peptides between enriched (PP) and
unenriched (FP) samples.

## Usage

``` r
checkRatioMat(ratioMat, minOverlap = 3)
```

## Arguments

- ratioMat:

  A numeric `matrix` representing the ratio of phosphoproteome data to
  full proteome data.

- minOverlap:

  A `numeric` specifying the minimum number of overlapping peptides
  required between samples. Default is 3.

## Value

A `character` vector of sample names that do not meet the overlap
criteria.
