# Plot Log Ratio of PP/FP (Phosphoproteome to Full Proteome) intensities

`plotLogRatio` generates a boxplot of the log2 ratio of intensities of
phosphoproteome to full proteome data from a `MultiAssayExperiment`
object.

## Usage

``` r
plotLogRatio(maeData, normalization = FALSE)
```

## Arguments

- maeData:

  A `MultiAssayExperiment` object containing phosphoproteome and full
  proteome data.

- normalization:

  A `logical` value indicating whether to perform normalization. Default
  is `FALSE`.

## Value

A `ggplot2` object representing the boxplot of the log2 ratios.

## Examples

``` r
# Load multiAssayExperiment object
data("dia_example")
# Call the function
plotLogRatio(dia_example, normalization = TRUE)

```
