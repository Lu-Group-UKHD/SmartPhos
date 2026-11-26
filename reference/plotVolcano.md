# Plot Volcano Plot for Differential Expression Analysis

`plotVolcano` generates a volcano plot to visualize differential
expression results.

## Usage

``` r
plotVolcano(tableDE, pFilter = 0.05, fcFilter = 0.5, usePadj = FALSE)
```

## Arguments

- tableDE:

  A `data frame` containing differential expression results with columns
  'ID', 'log2FC', 'pvalue', and 'Gene'.

- pFilter:

  A `numeric` value specifying the p-value threshold for significance.
  Default is 0.05.

- fcFilter:

  A `numeric` value specifying the log2 fold-change threshold for
  significance. Default is 0.5.

- usePadj:

  A `boolean` value specifying to use adjusted P-value or not. Default
  is FALSE.

## Value

A `ggplot2` object representing the volcano plot.

## Details

This function creates a volcano plot where differentially expressed
genes are categorized as 'Up', 'Down', or 'Not Sig' based on the
provided p-value and log2 fold-change thresholds. Points on the plot are
color-coded to indicate their expression status.

## Examples

``` r
library(SummarizedExperiment)
# Load multiAssayExperiment object
data("dda_example")
# Get SummarizedExperiment object
se <- dda_example[["Proteome"]]
colData(se) <- colData(dda_example)
# Preprocess the proteome assay
result <- preprocessProteome(se, normalize = TRUE)
# Call the function to perform differential expression analyis
de <- performDifferentialExp(se = result, assay = "Intensity",
method = "limma", reference = "1stCrtl", target = "EGF",
condition = "treatment")
# Plot the volcano plot from the result
plotVolcano(de$resDE)

```
