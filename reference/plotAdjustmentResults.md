# Plot Adjustment Results

`plotAdjustmentResults` generates plots to visualize the results of
phosphoproteome adjustment.

## Usage

``` r
plotAdjustmentResults(maeData, normalization = FALSE)
```

## Arguments

- maeData:

  A `MultiAssayExperiment` object containing phosphoproteome and
  proteome data.

- normalization:

  A `logical` value indicating whether normalization was performed.
  Default is `FALSE`.

## Value

A list containing:

- ratioTrendPlot:

  A `ggplot2` object showing the line plot of PP/FP ratios for features
  present in all samples.

- ratioBoxplot:

  A `ggplot2` object showing the box plot of PP/FP ratios before and
  after adjustment.

- ppBoxplot:

  A `ggplot2` object showing the box plot of phosphorylation intensities
  in PP samples before and after adjustment.

## Details

The function performs the following steps:

1.  Checks if the adjustment factor is present in the sample annotation.

2.  Calculates the ratio matrix before and after adjustment.

3.  Creates a trend line plot for features present in all samples.

4.  Creates box plots of the PP/FP ratios and phosphorylation
    intensities before and after adjustment.
