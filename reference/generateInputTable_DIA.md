# Generate Input Table for DIA Analysis

`generateInputTable_DIA` generates an input table for DIA analysis by
reading files from a specified folder.

## Usage

``` r
generateInputTable_DIA(rawFolder)
```

## Arguments

- rawFolder:

  A `character` string specifying the path to the folder containing the
  raw files.

## Value

A `data.frame` with columns fileName, searchType, and id that can be
used as input for further analysis.

## Details

The function performs the following steps:

- Reads the summary file containing experimental information.

- Generates unique experimental IDs based on sample type, treatment,
  timepoint, and replicate.

- Processes file paths for full proteome and phosphoproteome data.

- Creates a combined input table with file names, search types, and IDs.
