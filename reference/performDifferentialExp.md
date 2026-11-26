# Perform Differential Expression Analysis

`performDifferentialExp` performs differential expression analysis on a
given `SummarizedExperiment` object using either the `limma` or `ProDA`
method.

## Usage

``` r
performDifferentialExp(
  se,
  assay,
  method = c("limma", "ProDA"),
  condition = NULL,
  reference,
  target,
  refTime = NULL,
  targetTime = NULL,
  pairedTtest = FALSE
)
```

## Arguments

- se:

  A `SummarizedExperiment` object containing the data.

- assay:

  A `character` string specifying the assay to use for the analysis.

- method:

  A `character` string specifying the method to use for differential
  expression analysis ("limma" or "ProDA"). Default is "limma".

- condition:

  A `character` string specifying the condition column in colData(se).
  Default is `NULL`.

- reference:

  A `character` string or vector specifying the reference group.

- target:

  A `character` string or vector specifying the target group.

- refTime:

  A `character` string or vector specifying the reference time points.
  Default is `NULL`.

- targetTime:

  A `character` string or vector specifying the target time points.
  Default is `NULL`.

- pairedTtest:

  A `logical` value specifying to perform paired t-test or not. Default
  is `FALSE`.

## Value

A list containing:

- resDE:

  A `tibble` with the differential expression results.

- seSub:

  A `SummarizedExperiment` object subset to the samples used in the
  analysis.

## Details

This function is designed to facilitate differential expression analysis
on a `SummarizedExperiment` (SE) object. The function allows users to
specify various parameters to tailor the analysis to their specific
experimental setup.

The main steps of the function are as follows:

1\. Sample Selection: Based on the provided `condition`, `reference`,
and `target` arguments, the function identifies the relevant samples for
the analysis. If time points (`refTime` and `targetTime`) are provided,
it further refines the sample selection.

2\. Subsetting the SE Object: The SE object is subsetted to include only
the selected samples. A new column `comparison` is added to the
`colData`, indicating whether each sample belongs to the reference or
target group.

3\. Design Matrix Construction: The function constructs a design matrix
for the differential expression analysis. If the SE object contains a
`subjectID` column, this is included in the design to account for
repeated measures or paired samples.

4\. Differential Expression Analysis: Depending on the specified
`method`, the function performs the differential expression analysis
using either the `limma` or `ProDA` package: - `Limma`: The function
fits a linear model to the expression data and applies empirical Bayes
moderation to the standard errors. The results are then extracted and
formatted. - `ProDA`: The function fits a probabilistic dropout model to
the expression data and tests for differential expression. The results
are then extracted and formatted.

5\. Result Formatting: The differential expression results are merged
with the metadata from the SE object, and the resulting table is
formatted into a tibble. The table includes columns for log2 fold change
(log2FC), test statistic (stat), p-value (pvalue), adjusted p-value
(padj), and gene/feature ID (ID).

The function returns a `list` containing the formatted differential
expression results and the subsetted SE object. This allows users to
further explore or visualize the results as needed.

## Examples

``` r
library(SummarizedExperiment)
# Load multiAssayExperiment object
data("dda_example")
# Get SummarizedExperiment object
se <- dda_example[["Proteome"]]
colData(se) <- colData(dda_example)
# Preprocess the proteome assay
result <- preprocessProteome(se, normalize = TRUE)
# Call the function to perform differential expression analyis
performDifferentialExp(se = result, assay = "Intensity", method = "limma",
reference = "1stCrtl", target = "EGF", condition = "treatment")
#> $resDE
#> # A tibble: 83 × 8
#>    ID    log2FC  stat  pvalue  padj UniprotID            Gene     PeptideCounts
#>    <chr>  <dbl> <dbl>   <dbl> <dbl> <chr>                <chr>    <chr>        
#>  1 p99     2.02  4.75 0.00153 0.127 O00743               PPP6C    12           
#>  2 p95    -2.73 -4.03 0.00396 0.164 O00571               DDX3X    45           
#>  3 p48     2.53  3.73 0.00600 0.166 E9PRG8               C11orf98 4            
#>  4 p93     2.46  3.04 0.0165  0.276 O00560               SDCBP    10           
#>  5 p54    -2.23 -3.03 0.0166  0.276 O00151               PDLIM1   24           
#>  6 p38     2.83  2.43 0.0416  0.549 CON__P48668;P48668   KRT6C    46;46        
#>  7 p96     1.64  2.22 0.0577  0.549 O00622               CYR61    15           
#>  8 p100   -1.40 -2.21 0.0585  0.549 O00764               PDXK     19           
#>  9 p47    -2.03 -2.14 0.0655  0.549 E9PAV3;Q13765;Q9BZK3 NACA     9;9;2        
#> 10 p71    -1.55 -2.08 0.0720  0.549 O00299               CLIC1    15           
#> # ℹ 73 more rows
#> 
#> $seSub
#> class: SummarizedExperiment 
#> dim: 83 8 
#> metadata(0):
#> assays(1): Intensity
#> rownames(83): p2 p3 ... p99 p100
#> rowData names(3): UniprotID Gene PeptideCounts
#> colnames(8): Full_1stCrtl_0min_rep1 Full_1stCrtl_0min_rep2 ...
#>   Full_EGF_40min_rep1 Full_EGF_40min_rep2
#> colData names(7): sample sampleType ... batch comparison
#> 
```
