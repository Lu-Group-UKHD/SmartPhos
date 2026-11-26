# Read and Process a DIA Proteome Experiment

`readProteomeExperimentDIA` reads and processes DIA (Data-Independent
Acquisition) proteome data from multiple files and constructs a
`SummarizedExperiment` object.

## Usage

``` r
readProteomeExperimentDIA(fileTable, showProgressBar = FALSE)
```

## Arguments

- fileTable:

  A `data frame` containing metadata about the files to be read. Must
  contain columns searchType, fileName, id, and optionally outputID.

- showProgressBar:

  `Logical`, whether to show a progress bar during processing. Default
  is `FALSE`.

## Value

A `SummarizedExperiment` object containing the processed proteome data.
\#' @details The function performs the following steps:

- Filters the \`fileTable\` to include only rows where \`searchType\` is
  "proteome".

- For each file specified in \`fileTable\`, reads the data using
  \`data.table::fread\`.

- Removes rows where the \`PG.ProteinGroups\` column is NA or empty.

- Processes each sample in parallel using \`BiocParallel::bplapply\`,
  applying the \`readOneProteomDIA\` function to filter and clean the
  data for each sample.

- Combines the processed data from all files.

- Constructs a matrix of intensities with rows corresponding to proteins
  and columns corresponding to samples.

- Constructs a \`SummarizedExperiment\` object with the intensity matrix
  and protein annotations.

The `readOneProteomDIA` function is used to read and filter the data for
each individual sample, and it must be available in the environment.

## Examples

``` r
file <- system.file("extdata", "proteomeDIA_1.xls", package = "SmartPhos")
fileTable <- data.frame(searchType = "proteome", fileName = file,
id = c("sample1", "sample2"))
readProteomeExperimentDIA(fileTable)
#> class: SummarizedExperiment 
#> dim: 2 2 
#> metadata(0):
#> assays(1): Intensity
#> rownames(2): p1 p2
#> rowData names(2): UniprotID Gene
#> colnames(2): sample1 sample2
#> colData names(0):
```
