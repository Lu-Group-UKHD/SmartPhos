# Extract the Last Gene Symbol from a Semicolon-Separated List

`getOneSymbol` extracts the last gene symbol from a semicolon-separated
list of gene symbols.

## Usage

``` r
getOneSymbol(Gene)
```

## Arguments

- Gene:

  A `character` vector where each element is a semicolon-separated list
  of gene symbols.

## Value

A `character` vector containing the last gene symbol from each element
of the input vector.

## Details

This function processes a character vector where each element consists
of gene symbols separated by semicolons. It splits each element by
semicolons and extracts the last gene symbol from the resulting list.
The output is a character vector of these last gene symbols.
