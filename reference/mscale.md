# Scale and Center a Matrix

`mscale` scales and centers each row of a matrix, with options for using
mean or median, standard deviation or mean absolute deviation, and
censoring extreme values.

## Usage

``` r
mscale(x, center = TRUE, scale = TRUE, censor = NULL, useMad = FALSE)
```

## Arguments

- x:

  A `numeric` matrix where rows are features and columns are samples.

- center:

  `Logical`. If `TRUE`, the rows are centered by subtracting the mean or
  median. Default is `TRUE`.

- scale:

  `Logical`. If `TRUE`, the rows are scaled by dividing by the standard
  deviation or mean absolute deviation. Default is `TRUE`.

- censor:

  A `numeric` vector of length one or two for censoring the scaled
  values. If length one, values are censored symmetrically at positive
  and negative values. If length two, the first value is the lower limit
  and the second value is the upper limit. Default is `NULL`.

- useMad:

  `Logical`. If `TRUE`, the mean absolute deviation is used for scaling
  instead of the standard deviation. Default is `FALSE`.

## Value

A scaled and centered numeric `matrix` with the same dimensions as the
input `matrix` \`x\`.

## Details

The function allows for flexible scaling and centering of the rows of a
`matrix`:

- If both `center` and `scale` are `TRUE`, rows are centered and scaled.

- If only `center` is `TRUE`, rows are centered but not scaled.

- If only `scale` is `TRUE`, rows are scaled but not centered.

- If neither `center` nor `scale` is `TRUE`, the original `matrix` is
  returned.

The function can also censor extreme values, either symmetrically or
asymmetrically, based on the `censor` parameter.

## Examples

``` r
# Create a sample matrix (3 rows by 5 columns)
sample_matrix <- matrix(c(1:15), nrow = 3, byrow = TRUE)

# Scale and center the matrix using the default settings
mscale(sample_matrix, center = TRUE, scale = TRUE)
#>           [,1]       [,2] [,3]      [,4]     [,5]
#> [1,] -1.264911 -0.6324555    0 0.6324555 1.264911
#> [2,] -1.264911 -0.6324555    0 0.6324555 1.264911
#> [3,] -1.264911 -0.6324555    0 0.6324555 1.264911

# Only center the matrix without scaling
mscale(sample_matrix, center = TRUE, scale = FALSE)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]   -2   -1    0    1    2
#> [2,]   -2   -1    0    1    2
#> [3,]   -2   -1    0    1    2

# Only scale the matrix without centering
mscale(sample_matrix, center = FALSE, scale = TRUE)
#>           [,1]     [,2]     [,3]     [,4]     [,5]
#> [1,] 0.6324555 1.264911 1.897367 2.529822 3.162278
#> [2,] 3.7947332 4.427189 5.059644 5.692100 6.324555
#> [3,] 6.9570109 7.589466 8.221922 8.854377 9.486833
```
