# Normalize a Matrix Using Median or Mean

`medianNorm` normalizes the columns of a matrix by either the median or
the mean.

## Usage

``` r
medianNorm(x, method = "median")
```

## Arguments

- x:

  A `numeric` matrix to be normalized.

- method:

  A `character` string specifying the normalization method. Options are
  "median" or "mean". Default is "median".

## Value

A `numeric` `matrix` with normalized columns.

## Details

The function performs the following steps:

1.  If the `method` is "median", it calculates the median of each column
    and adjusts by the overall median of these medians.

2.  If the `method` is "mean", it calculates the mean of each column and
    adjusts by the overall mean of these means.

3.  It constructs a matrix of these adjusted values and subtracts it
    from the original matrix to normalize the columns.

## Examples

``` r
# Example usage:
x <- matrix(rnorm(20), nrow=5, ncol=4)
medianNorm(x, method = "median")
#>               [,1]       [,2]       [,3]       [,4]
#> [1,]  0.7509384756  0.0444788 -0.4219140 -0.5040665
#> [2,]  0.0004157175  0.2717034  1.2377762  0.1361387
#> [3,]  0.1116497372  0.1671684  0.1116497  0.1116497
#> [4,]  1.9740924522  0.1116497 -0.3420521  1.5082484
#> [5,] -0.0366068151 -1.2361092  0.9433930 -0.5291011
```
