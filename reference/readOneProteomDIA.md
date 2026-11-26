# Read and Process a Single DIA Proteomics Sample

`readOneProteomDIA` reads and processes data from a single DIA
proteomics sample, applying filtering and data transformation steps.

## Usage

``` r
readOneProteomDIA(inputTab, sampleName)
```

## Arguments

- inputTab:

  A `data.table` or `data.frame` containing the input data.

- sampleName:

  A `character` string specifying the sample name.

## Value

A `data.table` containing the processed data with columns for intensity,
UniProt ID, and gene name.

## Details

This function processes DIA proteomics data for a single sample by
filtering out rows with non-quantitative data, converting character
values to numeric, and renaming columns for consistency. It also ensures
that each protein group has a unique identifier.
