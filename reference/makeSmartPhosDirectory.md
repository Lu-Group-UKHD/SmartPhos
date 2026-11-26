# Create SmartPhos Directory Structure

`makeSmartPhosDirectory` creates a directory for the SmartPhos shiny
app, and copies the necessary Shiny app files into the newly created
directory.

## Usage

``` r
makeSmartPhosDirectory(path)
```

## Arguments

- path:

  A `character` string specifying the directory path where the SmartPhos
  folder should be created.

## Value

None (invisible NULL). The function creates the necessary directories
and copies files.

## Details

The function first creates the main directory at the specified path and
a subdirectory named \`"save"\` for storing `MultiAssayExperiment`
object. It then locates the Shiny application files from the SmartPhos
package and copies them into the new directory.

## Examples

``` r
makeSmartPhosDirectory("shinyApp")
#> Successful!!
```
