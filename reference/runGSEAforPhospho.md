# Run GSEA for Phosphorylation Data

`runGSEAforPhospho` performs Gene Set Enrichment Analysis (GSEA) for
phosphorylation data.

## Usage

``` r
runGSEAforPhospho(
  geneStat,
  ptmSetDb,
  nPerm,
  weight = 1,
  correl.type = c("rank", "symm.rank", "z.score"),
  statistic = c("Kolmogorov-Smirnov", "area.under.RES"),
  min.overlap = 5
)
```

## Arguments

- geneStat:

  A `data frame` containing gene statistics, with gene names as row
  names and a column named 'stat' for the statistics.

- ptmSetDb:

  A `data frame` of post-translational modification (PTM) signature
  sets.

- nPerm:

  A `numeric` value specifying the number of permutations for the null
  distribution.

- weight:

  A `numeric` value for the weight parameter in the GSEA algorithm. If
  weight == 0 then the test statistics do not matter. Default is 1.

- correl.type:

  A `character` string specifying the correlation type. Options are
  "rank", "symm.rank", and "z.score". Default is "rank".

- statistic:

  A `character` string specifying the statistic to be used. Options are
  "Kolmogorov-Smirnov" and "area.under.RES". Default is
  "Kolmogorov-Smirnov".

- min.overlap:

  A `numeric` specifying the minimum overlap required between gene sets
  and the input data. Default is 5.

## Value

A `tibble` with enrichment scores and associated statistics for each PTM
set.

## Details

This function runs GSEA on phosphorylation data to identify enriched PTM
sets. It calculates enrichment scores and p-values for each set,
normalizes the scores, and adjusts p-values for multiple testing.
