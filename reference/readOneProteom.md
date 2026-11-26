# Read and Process One Proteomics Sample

`readOneProteom` reads and processes proteomics data for a single
sample, applying filters for peptide count and optionally using LFQ
quantification. It returns a `data.table` with useful columns and unique
identifiers.

## Usage

``` r
readOneProteom(inputTab, sampleName, pepNumCut = 1, ifLFQ = TRUE)
```

## Arguments

- inputTab:

  A `data.table` or `data.frame` containing the input data for the
  proteomics sample.

- sampleName:

  A `character` string specifying the name of the sample to be
  processed.

- pepNumCut:

  A `numeric` value specifying the minimum number of peptides required
  for a protein to be included. Default is 1.

- ifLFQ:

  A `logical` value indicating whether to use LFQ quantification.
  Default is `TRUE`.

## Value

A `data.table` with the processed proteomics data, including columns for
intensity, Uniprot ID, peptide counts, and gene names.

## Details

This function processes proteomics data for a single sample by filtering
based on the number of peptides and optionally using LFQ quantification.
It ensures that unique identifiers are created for each protein, and
removes rows with missing or zero quantification values.
