# Plot Kinase score for Differential Expression data

\`plotKinaseDE\` generates a `bar plot` of the top kinases associated
with the differentially expressed genes based on their scores.

## Usage

``` r
plotKinaseDE(scoreTab, nTop = 10, pCut = 0.05)
```

## Arguments

- scoreTab:

  A `data frame` containing kinase scores with columns source, score,
  and p_value.

- nTop:

  A `numeric` value specifying the number of top kinases to plot for
  each direction. Default is 10.

- pCut:

  A `numeric` value specifying the p-value cutoff for significance.
  Default is 0.05.

## Value

A `ggplot2` object representing the bar plot of kinase score.

## Details

The function performs the following steps:

1.  Adds a column for significance based on the p-value cutoff.

2.  Adds a column for the sign of the score.

3.  Filters out kinases with a score of 0.

4.  Selects the top `nTop` kinases by absolute score for each sign of
    the score.

5.  Creates a bar plot with the selected kinases.

## Examples

``` r
# Example usage:
scoreTab <- data.frame(
 source = c("Kinase1", "Kinase2", "Kinase3", "Kinase4"),
 score = c(2.3, -1.5, 0, 3.1),
 p_value = c(0.01, 0.2, 0.05, 0.03)
)
plotKinaseDE(scoreTab, nTop = 3, pCut = 0.05)

```
