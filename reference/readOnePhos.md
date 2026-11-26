# Read and Filter Phosphorylation Data for a Specific Sample

`readOnePhos` reads phosphorylation data from an input table, filters it
based on localization probability, score difference, and intensity, and
returns the filtered data for a specific sample.

## Usage

``` r
readOnePhos(
  inputTab,
  sampleName,
  localProbCut = 0.75,
  scoreDiffCut = 5,
  multiMap
)
```

## Arguments

- inputTab:

  A `data.table` or `data.frame` containing phosphorylation data with
  columns for localization probability, score difference, and intensity
  for various samples.

- sampleName:

  A `character` string specifying the sample name to filter data for.

- localProbCut:

  A `numeric` value specifying the cutoff for localization probability.
  Default is 0.75.

- scoreDiffCut:

  A `numeric` value specifying the cutoff for score difference. Default
  is 5.

- multiMap:

  A `logical` value indicating whether to allow multiple mapping (not
  used in this function but could be relevant for further extensions).

## Value

A `data.frame` containing the filtered phosphorylation data for the
specified sample, with columns for intensity, Uniprot ID, gene name,
position within proteins, amino acid residue, and sequence window.

## Details

The function filters the input phosphorylation data based on three
criteria: localization probability, score difference, and intensity.
Only rows that meet or exceed the specified cutoffs for these criteria
and have non-zero intensity are retained. The filtered data is then
returned with a unique identifier for each row.
