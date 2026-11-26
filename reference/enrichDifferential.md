# Perform Enrichment analysis on differentially expressed genes or phospho-sites

`enrichDifferential` performs enrichment analysis on differentially
expressed genes and phospho-sites for either pathway or phospho-specific
enrichment, depending on the input parameters. It supports multiple
statistical methods such as PAGE and GSEA for pathway enrichment and a
Kolmogorov-Smirnov approach for phospho-enrichment.

## Usage

``` r
enrichDifferential(
  dea,
  type = c("Pathway enrichment", "Phospho-signature enrichment"),
  gsaMethod = c("PAGE", "GSEA"),
  geneSet,
  ptmSet,
  statType = c("stat", "log2FC"),
  nPerm = 100,
  sigLevel = 0.05,
  ifFDR = FALSE
)
```

## Arguments

- dea:

  A `data frame` containing the differential expression analysis
  results. It should include columns like \`pvalue\`, \`Gene\` (or
  \`site\`), \`stat\`, and \`log2FC\`.

- type:

  A `character` string indicating the type of enrichment. Options are
  \`"Pathway enrichment"\` or \`"Phospho-signature enrichment"\`.

- gsaMethod:

  A `character` string specifying the gene set analysis method for
  pathway enrichment. Options are \`"PAGE"\` or \`"GSEA"\`.

- geneSet:

  A gene set collection to use for pathway enrichment.

- ptmSet:

  A post-translational modification (PTM) set database for
  phospho-enrichment analysis.

- statType:

  A `character` string specifying the statistic type to use. Options are
  \`"stat"\` or \`"log2FC"\`.

- nPerm:

  A `numeric` specifying the number of permutations for GSEA. Default is
  100.

- sigLevel:

  A `numeric` value representing the significance threshold for
  filtering results. Ddefault is 0.05.

- ifFDR:

  A `logical` value indicating whether to filter results using
  FDR-adjusted p-values. Default is \`FALSE\`.

## Value

A data frame containing the results of the enrichment analysis,
including columns such as the gene set name, statistical significance,
and adjusted p-values.

## Details

The \`enrichDifferential\` function performs either pathway enrichment
or phospho-enrichment analysis based on the \`type\` parameter. For
pathway enrichment, it uses either the PAGE or GSEA method with a
provided gene set collection. For phospho-enrichment, it uses a
Kolmogorov-Smirnov test with a PTM set database. Results can be filtered
by significance level and optionally adjusted for FDR.

## Examples

``` r
library(SummarizedExperiment)
library(piano)
# Load multiAssayExperiment object
data("dia_example")
# Get SummarizedExperiment object
se <- dia_example[["Phosphoproteome"]]
colData(se) <- colData(dia_example)
# Preprocess the proteome assay
result <- preprocessPhos(se, normalize = TRUE)
# Call the function to perform differential expression analyis
dea <- performDifferentialExp(se = result, assay = "Intensity",
method = "limma", reference = "1stCrtl", target = "EGF",
condition = "treatment")
# Load the gene set
genesetPath <- appDir <- system.file("shiny-app/geneset",
package = "SmartPhos")
inGMT <- loadGSC(paste0(genesetPath,"/Cancer_Hallmark.gmt"),type="gmt")
# Call the function
resTab <- enrichDifferential(dea = dea$resDE, type = "Pathway enrichment",
gsaMethod = "PAGE", geneSet = inGMT, statType = "stat", nPerm = 200,
sigLevel = 0.05, ifFDR = FALSE)
#> Running gene set analysis:
#> Checking arguments...done!
#> Final gene/gene-set association: 20 genes and 22 gene sets
#>   Details:
#>   Calculating gene set statistics from 20 out of 65 gene-level statistics
#>   Removed 4366 genes from GSC due to lack of matching gene statistics
#>   Removed 28 gene sets containing no genes after gene removal
#>   Removed additionally 0 gene sets not matching the size limits
#>   Loaded additional information for 22 gene sets
#> Calculating gene set statistics...done!
#> Calculating gene set significance...done!
#> Adjusting for multiple testing...done!
```
