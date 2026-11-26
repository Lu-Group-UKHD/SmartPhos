# Plot PCA

`plotPCA` generates a PCA plot using the results from a PCA analysis and
a `SummarizedExperiment` object. The points on the plot can be colored
and shaped based on metadata.

## Usage

``` r
plotPCA(pca, se, xaxis = "PC1", yaxis = "PC2", color = "none", shape = "none")
```

## Arguments

- pca:

  A PCA result object, typically obtained from `prcomp`.

- se:

  A `SummarizedExperiment` object containing the metadata.

- xaxis:

  A `character` string specifying which principal component to use for
  the x-axis. Default is "PC1".

- yaxis:

  A `character` string specifying which principal component to use for
  the y-axis. Default is "PC2".

- color:

  A `character` string specifying the metadata column to use for
  coloring the points. Default is "none".

- shape:

  A `character` string specifying the metadata column to use for shaping
  the points. Default is "none".

## Value

A `ggplot2` object showing the PCA plot.

## Details

This function creates a PCA plot using the scores from a PCA result
object and metadata from a `SummarizedExperiment` object. The x-axis and
y-axis can be customized to display different principal components, and
the points can be optionally colored and shaped based on specified
metadata columns.

## Examples

``` r
# Load multiAssayExperiment object
data("dia_example")
# Get SummarizedExperiment object
se <- dia_example[["Phosphoproteome"]]
SummarizedExperiment::colData(se) <- SummarizedExperiment::colData(
dia_example)
# Generate the imputed assay
result <- preprocessPhos(seData = se, normalize = TRUE, impute = "QRILC")
#> Imputing along margin 2 (samples/columns).
# Perform PCA
pcaResult <- stats::prcomp(t(
SummarizedExperiment::assays(result)[["imputed"]]),
center = TRUE, scale. = TRUE)
# Plot PCA results
plotPCA(pca = pcaResult, se = result, color = "treatment")
#> Error: unable to find an inherited method for function ‘plotPCA’ for signature ‘object = "missing"’
```
