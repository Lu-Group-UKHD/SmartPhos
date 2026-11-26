# Launch the SmartPhos Shiny Application

`runSmartPhos` launches the SmartPhos `Shiny` application, which
provides an interactive interface for analyzing phosphoproteomic data.

## Usage

``` r
runSmartPhos()
```

## Value

The function does not return a value; it starts the Shiny application
for `SmartPhos`.

## Details

The `runSmartPhos` function locates the `Shiny` app directory within the
`SmartPhos` package and launches the application. If the app directory
cannot be found, the function will stop and prompt the user to
re-install the `SmartPhos` package.

## Examples

``` r
# To run the SmartPhos Shiny application, simply call:
# runSmartPhos()
```
