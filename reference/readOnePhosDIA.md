# Read Phosphorylation Data for One Sample from DIA

`readOnePhosDIA` reads and processes phosphorylation data for a single
sample from a DIA experiment, applying filters for localization
probability and removing duplicates if specified.

## Usage

``` r
readOnePhosDIA(inputTab, sampleName, localProbCut = 0.75, removeDup = FALSE)
```

## Arguments

- inputTab:

  A `data.table` or `data.frame` containing phosphorylation data.

- sampleName:

  A `character` string specifying the sample name.

- localProbCut:

  A `numeric` value specifying the cutoff for localization probability.
  Default is 0.75.

- removeDup:

  A `logical` value indicating whether to remove duplicate entries based
  on UniprotID and intensity. Default is `FALSE`.

## Value

A `data.table` containing the processed phosphorylation data for the
specified sample.

## Details

This function processes phosphorylation data for a single sample by
filtering based on localization probability and non-zero intensity. It
handles multiplicity by summarizing intensities and optionally removes
duplicates. The resulting data is returned as a data.table with unique
identifiers.
