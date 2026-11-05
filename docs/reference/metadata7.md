# Get or set experiment-wide metadata

Get or set experiment-wide metadata

## Usage

``` r
metadata7(x)

metadata7(x) <- value
```

## Arguments

- x:

  A `SummarizedExperiment7`.

- value:

  Replacement named list.

## Examples

``` r
se <- SummarizedExperiment7(assays = list(counts = matrix(1:6, 2, 3)))
metadata7(se) <- list(study = "demo")
metadata7(se)
#> $study
#> [1] "demo"
#> 
```
