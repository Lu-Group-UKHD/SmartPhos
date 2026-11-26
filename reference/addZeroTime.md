# Add Zero Timepoint Data to Treatment Subset

`addZeroTime` adds a zero timepoint to a specific treatment's data
subset.

## Usage

``` r
addZeroTime(data, condition, treat, zeroTreat, timeRange)
```

## Arguments

- data:

  A `SummarizedExperiment` object containing the experimental data.

- condition:

  `Character` string corresponds to one of the columns from the colData
  of SE object.

- treat:

  `Character` string specifying the treatment to which zero timepoint
  should be added.

- zeroTreat:

  `Character` string specifying the treatment representing the zero
  timepoint.

- timeRange:

  `Character` vector specifying the timepoints to include for the
  treatment.

## Value

A `SummarizedExperiment` object with the zero timepoint added to the
specified treatment's data.

## Details

The function performs the following steps:

1.  Subsets the data for the specified treatment and time range.

2.  Subsets the data for the zero timepoint of the specified zero
    treatment.

3.  Combines the assays from the treatment and zero timepoint subsets.

4.  Updates the column data to reflect the combined treatment.

5.  Returns a `SummarizedExperiment` object with the combined data.

## Examples

``` r
library(SummarizedExperiment)
#> Loading required package: MatrixGenerics
#> Loading required package: matrixStats
#> 
#> Attaching package: ‘MatrixGenerics’
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
#>     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
#>     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
#>     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
#>     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
#>     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
#>     colWeightedMeans, colWeightedMedians, colWeightedSds,
#>     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
#>     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
#>     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
#>     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
#>     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
#>     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
#>     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
#>     rowWeightedSds, rowWeightedVars
#> Loading required package: GenomicRanges
#> Loading required package: stats4
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following object is masked from ‘package:SmartPhos’:
#> 
#>     plotPCA
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: S4Vectors
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> Loading required package: IRanges
#> Loading required package: Seqinfo
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> 
#> Attaching package: ‘Biobase’
#> The following object is masked from ‘package:MatrixGenerics’:
#> 
#>     rowMedians
#> The following objects are masked from ‘package:matrixStats’:
#> 
#>     anyMissing, rowMedians
# Load multiAssayExperiment object
data("dia_example")
# Get SummarizedExperiment object
se <- dia_example[["Phosphoproteome"]]
colData(se) <- colData(dia_example)
# Call the function
addZeroTime(se, condition = "treatment", treat = "EGF",
zeroTreat = "1stCrtl", timeRange = c("20min","40min", "6h"))
#> class: SummarizedExperiment 
#> dim: 500 6 
#> metadata(0):
#> assays(1): intensity
#> rownames(500): s1 s2 ... s499 s500
#> rowData names(7): UniprotID Gene ... Sequence site
#> colnames(6): FullProteome_EGF_0min_rep2 FullProteome_EGF_0min_rep3 ...
#>   FullProteome_EGF_0min_rep1 Phospho_EGF_0min_rep1
#> colData names(6): sample treatment ... sampleType sampleName
```
