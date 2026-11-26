# Plot Heatmap of Intensity assay

`plotHeatmap` generates a heatmap for intensity assay for different
conditions, including top variants, differentially expressed genes, and
selected time series clusters.

## Usage

``` r
plotHeatmap(
  type = c("Top variant", "Differentially expressed", "Selected time series cluster"),
  se,
  data = NULL,
  top = 100,
  cutCol = 1,
  cutRow = 1,
  clustCol = TRUE,
  clustRow = TRUE,
  annotationCol = NULL,
  title = NULL
)
```

## Arguments

- type:

  A `character` string specifying the type of heatmap to plot. Options
  are "Top variant", "Differentially expressed", and "Selected time
  series cluster".

- se:

  A `SummarizedExperiment` object containing the imputed intensity
  assay.

- data:

  An optional `data frame` containing additional data for
  "Differentially expressed" and "Selected time series cluster" types.
  Default is `NULL`.

- top:

  A `numeric` value specifying the number of top variants to plot.
  Default is 100.

- cutCol:

  A `numeric` value specifying the number of clusters for columns.
  Default is 1.

- cutRow:

  A `numeric` value specifying the number of clusters for rows. Default
  is 1.

- clustCol:

  A `logical` value indicating whether to cluster columns. Default is
  `TRUE`.

- clustRow:

  A `logical` value indicating whether to cluster rows. Default is
  `TRUE`.

- annotationCol:

  A `character` vector specifying the columns in the metadata to use for
  annotation. Default is `NULL`.

- title:

  A `character` string specifying the title of the heatmap. Default is
  `NULL`.

## Value

A `pheatmap` object showing the heatmap of Intensity data.

## Details

This function creates a heatmap using the Intensity assay from a
`SummarizedExperiment` object. The heatmap can show the top variants
based on standard deviation, differentially expressed genes, or selected
time series clusters. Row normalization is performed, and the heatmap
can include annotations based on specified metadata columns.

## Examples

``` r
library(SummarizedExperiment)
# Load multiAssayExperiment object
data("dia_example")
# Get SummarizedExperiment object
se <- dia_example[["Phosphoproteome"]]
colData(se) <- colData(dia_example)
# Generate the imputed assay
result <- preprocessPhos(seData = se, normalize = TRUE, impute = "QRILC")
#> Imputing along margin 2 (samples/columns).
# Plot heatmap for top variant
plotHeatmap(type = "Top variant", top = 10, se = result, cutCol = 2)

```
