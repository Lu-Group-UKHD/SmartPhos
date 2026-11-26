# Generate Input Table for Proteomic and Phosphoproteomic Analysis

`generateInputTable` generates an input table for proteomic and
phosphoproteomic analysis by reading files from a specified folder.

## Usage

``` r
generateInputTable(rawFolder, batchAsFolder = FALSE)
```

## Arguments

- rawFolder:

  A `character` string specifying the path to the folder containing the
  raw files.

- batchAsFolder:

  A `logical` value indicating whether to treat subdirectories as
  separate batches. Default is `FALSE`.

## Value

A `data.frame` with columns fileName, sample, searchType, batch, and id
that can be used as input for further analysis.

## Details

The function performs the following steps:

- Optionally treats subdirectories as separate batches.

- Reads the summary file containing experimental information.

- Generates unique experimental IDs based on batch and sample names.

- Processes file paths for total proteome / full proteome and
  phosphoproteome data.

- Creates a combined input table with file names, sample names, search
  types, batches, and IDs.
